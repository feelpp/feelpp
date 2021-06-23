//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/feelcrb/crbplugin.hpp>
#include "electric.hpp"

namespace Feel {

Electric::Electric()
    : super_type( "electric" )
{}

Electric::Electric( mesh_ptrtype mesh )
    : super_type( "electric" )
{
    this->M_mesh = mesh;
}

int Electric::Qa()
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int matSize = materials.materialWithPhysic("electric").size();
    return matSize+dirSize;
}

int Electric::mQA( int q )
{
    return 1;
}

int Electric::Nl()
{
    return 2;
}

int Electric::Ql( int l)
{
    auto bc = M_modelProps->boundaryConditions();
    int sourceSize = bc["potential"]["SourceTerm"].size();
    int dirSize = bc["potential"]["Dirichlet"].size();
    switch(l) {
    case 0:
        return sourceSize+dirSize;
    case 1:
        return 1;
    default:
        return 0;
    }
}

int Electric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    case 1:
        return mIntensity(q);
    default:
        return 0;
    }
}

int Electric::mCompliantQ(int q )
{
    return 1;
}

int Electric::mIntensity(int q )
{
    return 1;
}

void Electric::resizeQm( bool resizeMatrix )
{
    if( resizeMatrix )
        M_Aqm.resize( Qa());
    M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        if( resizeMatrix )
            M_Aqm[q].resize(mQA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mQA(q));
    }

    if( resizeMatrix )
        M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        if( resizeMatrix )
            M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            if( resizeMatrix )
                M_Fqm[l][q].resize(mLQF(l, q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mLQF(l, q) );
        }
    }

    if( resizeMatrix )
    {
        M_InitialGuess.resize(1);
        M_InitialGuess[0].resize(1);
        M_InitialGuess[0][0] = Xh->elementPtr();
    }
}

Electric::parameter_type
Electric::paramFromProperties() const
{
    auto mu = Dmu->element();
    int i = 0;
    auto parameters = M_modelProps->parameters();
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu(i++) = parameterPair.second.value();
        }
    }
    return mu;
}

void Electric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    std::string propertyPath = Environment::expand( soption("thermoelectric.filename"));
    M_modelProps = std::make_shared<prop_type>(propertyPath);
    M_modelProps->enableBoundaryConditions2();
    this->addModelFile("property-file", propertyPath);

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu->setDimension(nbCrbParameters);
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    std::vector<std::string> range;
    auto materials = M_modelProps->materials();
    auto electroMat = materials.materialWithPhysic("electric");
    for( auto const& mat : electroMat )
        range.push_back(mat.first);

    auto domain = markedelements(M_mesh, range);
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=domain ) );

    Feel::cout << "Potential nDof  : " << Xh->nDof() << std::endl;

    if( !M_V )
        M_V = element_ptrtype( new element_type( Xh ) );

    Jh = J_space_type::New( _mesh=M_mesh, _range=domain );

    this->resizeQm();
    this->decomposition();
}

void Electric::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    Feel::cout << "setupSpecificityModel" << std::endl;

    std::string propertyPath;
    if( this->hasModelFile("property-file") )
        propertyPath = this->additionalModelFiles().find("property-file")->second;
    else
        Feel::cerr << "Warning!! the database does not contain the property file! Expect bugs!"
                   << std::endl;
    M_modelProps = std::make_shared<prop_type>(propertyPath);

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu = parameterspace_type::New( nbCrbParameters, Environment::worldComm() );

    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    this->resizeQm(false);
}

void Electric::decomposition()
{
    auto V = Xh->element();
    auto phiV = Xh->element();

    auto gamma = doption("thermoelectric.gamma");

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("electric");

    /************** Left hand side **************/
    // electro
    int idx = 0;
    for( auto const& mat : materials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(V),grad(phiV)) );
        M_Aqm[idx++][0] = a0.matrixPtr();
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto aVD = form2(_test=Xh, _trial=Xh);
        aVD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         gamma/hFace()*inner(idt(V),id(phiV))
                         -inner(gradt(V)*N(),id(phiV))
                         -inner(grad(phiV)*N(),idt(V)) );
        M_Aqm[idx++][0] = aVD.matrixPtr();
    }

    /************** Right hand side **************/
    idx = 0;
    for( auto const& exAtM : bc["potential"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression2());
        auto fVST = form1(_test=Xh);
        fVST = integrate( markedelements(M_mesh, exAtM.marker() ),
                          e*id(phiV) );
        M_Fqm[0][idx++][0] = fVST.vectorPtr();
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto fVD = form1(_test=Xh);
        fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         gamma/hFace()*id(phiV) -  grad(phiV)*N() );
        M_Fqm[0][idx++][0] = fVD.vectorPtr();
    }

    // Intensity output
    auto fI = form1(_test=Xh);
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( exAtM.expression() != "0" )
            fI += integrate( markedfaces(M_mesh, exAtM.marker() ),
                             -grad(phiV)*N() );
    }
    M_Fqm[1][0][0] = fI.vectorPtr();

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(V),grad(phiV)) );
    M_energy_matrix = m.matrixPtr();
}

Electric::beta_vector_type
Electric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

Electric::beta_type
Electric::computeBetaQm( element_type const& T, parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

Electric::beta_type
Electric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

Electric::beta_type
Electric::computeBetaQm( parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void Electric::fillBetaQm( parameter_type const& mu )
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("electric");

    // Left Hand Side
    int idx = 0;
    for( auto const& mat : materials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("sigmaKey"));
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        M_betaAqm[idx++][0] = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));

    // Right Hand Side
    idx = 0;
    for( auto const& exAtM : bc["potential"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[0][idx++][0] = e.evaluate();
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));

        M_betaFqm[0][idx++][0] = sigma*e.evaluate();
    }

    // Intensity Output
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( exAtM.expression() != "0" )
            M_betaFqm[1][0][0] = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
    }
}

Electric::beta_vector_type
Electric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

Electric::affine_decomposition_type
Electric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<Electric::sparse_matrix_ptrtype> >
Electric::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<Electric::element_ptrtype> >
Electric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

Electric::element_type
Electric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto V = Xh->element();
    auto phiV = Xh->element();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("electric");
    auto gamma = doption("thermoelectric.gamma");
    auto current = mu.parameterNamed("current");

    /***************************** Electro *****************************/
    tic();
    auto a = form2(_test=Xh, _trial=Xh);
    // V
    for( auto const& mat : materials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        a += integrate( markedelements(M_mesh, mat.first),
                        sigma*inner(gradt(V),grad(phiV)) );
    }

    // V Dirichlet condition
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
        a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                        sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                               -inner(gradt(V)*N(),id(phiV))
                               -inner(grad(phiV)*N(),idt(V)) ) );
    }

    auto f = form1(_test=Xh);
    for( auto const& exAtM : bc["potential"]["SourceTerm"] )
    {
        auto e1 = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e1.expression().hasSymbol(param.first) )
                e1.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto e2 = expr(exAtM.expression2());

        f += integrate( markedelements(M_mesh, exAtM.marker() ),
                        e2*e1.evaluate()*id(phiV) );
    }
    // V Dirichlet condition
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));

        f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                        sigma*e.evaluate()*(gamma/hFace()*id(phiV) -  grad(phiV)*N()) );
    }
    a.solve( _solution=V, _rhs=f, _name="electro" );

    if( boption("thermoelectric.export-FE") )
    {
        auto e = exporter(M_mesh);
        auto j = Jh->element();
        for( auto const& mat : materials )
        {
            auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
            j += vf::project( _range=markedelements(M_mesh, mat.first), _space=Jh,
                              _expr=sigma*inner(gradv(V),gradv(V)) );
        }
        e->add( "joule", j );
        e->add("V", V);
        e->save();
    }
    toc("mono");

    return V;
}

double Electric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index == 0 || output_index == 1 )
    {
        for ( int q = 0; q < Ql(output_index); q++ )
        {
            element_ptrtype eltF( new element_type( Xh ) );
            *eltF = *M_Fqm[output_index][q][0];
            output += M_betaFqm[output_index][q][0]*dot( *eltF, u );
            //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

int Electric::mMaxSigma()
{
    return 1;
}

Electric::q_sigma_element_type Electric::eimSigmaQ(int m)
{
    q_sigma_element_type q = Xh->element();
    q.on( _range=elements(M_mesh), _expr=cst(1.) );
    return q;
}

Electric::vectorN_type Electric::eimSigmaBeta( parameter_type const& mu )
{
    vectorN_type beta(1);
    beta(0) = mu.parameterNamed( "sigma" );
    return beta;
}

void Electric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu )
{
    auto V = this->solve(mu);
    auto sigma = mu.parameterNamed("sigma");
    auto Vh = j.functionSpace();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
}

void Electric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu, element_type& V )
{
    auto sigma = mu.parameterNamed("sigma");
    auto Vh = j.functionSpace();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
}

FEELPP_CRB_PLUGIN( Electric, electric)
}
