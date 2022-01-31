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

#include <feel/feelmor/crbplugin.hpp>
#include "thermoelectric-simple.hpp"

namespace Feel {

ThermoElectric::ThermoElectric()
    : super_type( "thermoelectric-simple" )
{}

ThermoElectric::ThermoElectric( mesh_ptrtype mesh )
    : super_type( "thermoelectric-simple" )
{
    this->M_mesh = mesh;
}

int ThermoElectric::Qa()
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    return 2*materials.size() + dirSize + robSize;
}

int ThermoElectric::mQA( int q )
{
    return 1;
}

int ThermoElectric::Nl()
{
    return 3;
}

int ThermoElectric::Ql( int l)
{
    auto bc = M_modelProps->boundaryConditions();
    int sourceSize = bc["temperature"]["SourceTerm"].size();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    switch(l) {
    case 0:
        return sourceSize + dirSize + robSize;
    case 1:
        return 1;
    case 2:
        return 1;
    default:
        return 0;
    }
}

int ThermoElectric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    case 1:
        return mAverageTempQ(q);
    case 2:
        return mIntensityQ(q);
    default:
        return 0;
    }
}

int ThermoElectric::mCompliantQ(int q )
{
    return 1;
}

int ThermoElectric::mAverageTempQ(int q )
{
    return 1;
}

int ThermoElectric::mIntensityQ(int q )
{
    return 1;
}

void ThermoElectric::resizeQm( bool resizeMatrix )
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

ThermoElectric::parameter_type
ThermoElectric::paramFromProperties() const
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

void ThermoElectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    std::string propertyPath = Environment::expand( soption("thermoelectric.filename"));
    M_modelProps = std::make_shared<prop_type>(propertyPath);
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
    auto thermoMat = materials.materialWithPhysic("thermoelectric");
    for( auto const& mat : thermoMat )
        range.push_back(mat.first);

    auto domain = markedelements(M_mesh, range);
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=domain ) );

    Feel::cout << "Potential nDof  : " << Xh->template functionSpace<0>()->nDof() << std::endl
               << "Temperature nDof: " << Xh->template functionSpace<1>()->nDof() << std::endl;

    if( !M_VT )
        M_VT = element_ptrtype( new element_type( Xh ) );

    M_Jh = J_space_type::New( _mesh=M_mesh, _range=domain );

    this->resizeQm();
    this->decomposition();
}

void ThermoElectric::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
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

void ThermoElectric::decomposition()
{
    auto VT = Xh->element();
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();
    auto phi = Xh->element();
    auto phiV = phi.template element<0>();
    auto phiT = phi.template element<1>();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    auto gamma = doption("thermoelectric.gamma");

    /************** Left hand side **************/
    int idx = 0;
    // electro
    for( auto const& mat : materials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(V),grad(phiV)) );
        M_Aqm[idx++][0] = a0.matrixPtr();
    }
    // thermo
    for( auto const& mat : materials )
    {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(T),grad(phiT)) );
        M_Aqm[idx++][0] = a1.matrixPtr();
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
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto aTR = form2(_test=Xh, _trial=Xh);
        aTR += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          inner(idt(T),id(phiT)) );
        M_Aqm[idx++][0] = aTR.matrixPtr();
    }

    /************** Right hand side **************/
    int output = 0;
    idx = 0;
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression2());
        auto fTST = form1(_test=Xh);
        fTST = integrate( markedelements(M_mesh, exAtM.marker() ),
                          e*id(phiT) );
        M_Fqm[output][idx++][0] = fTST.vectorPtr();
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto fVD = form1(_test=Xh);
        fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         gamma/hFace()*id(phiV) -  grad(phiV)*N() );
        M_Fqm[output][idx++][0] = fVD.vectorPtr();
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto fTR = form1(_test=Xh);
        fTR = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         id(phiT) );
        M_Fqm[output][idx++][0] = fTR.vectorPtr();
    }
    output++;

    // Average temperature output
    auto fAvgT = form1(_test=Xh);
    for( auto const& mat : materials )
    {
        double area = integrate( markedelements(M_mesh, mat.first), cst(1.0) ).evaluate()(0,0) ;
        fAvgT += integrate( markedelements(M_mesh, mat.first ),
                            id(phiT)/cst(area) );
    }
    M_Fqm[output++][0][0] = fAvgT.vectorPtr();

    // Intensity output
    auto fI = form1(_test=Xh);
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( exAtM.expression() != "0" )
            fI += integrate( markedfaces(M_mesh, exAtM.marker() ),
                             -grad(phiV)*N() );
    }
    M_Fqm[output++][0][0] = fI.vectorPtr();
    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(V),grad(phiV)) + inner(gradt(T),grad(phiT)) );
    M_energy_matrix = m.matrixPtr();
}

ThermoElectric::beta_vector_type
ThermoElectric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( element_type const& T, parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void ThermoElectric::fillBetaQm( parameter_type const& mu )
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");

    // Left Hand Side
    int idx = 0;
    for( auto const& mat : materials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("sigmaKey"));
    for( auto const& mat : materials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("kKey"));
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        M_betaAqm[idx++][0] = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaAqm[idx++][0] = e.evaluate();
    }

    // Right Hand Side
    int output = 0;
    idx = 0;
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate();
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));

        M_betaFqm[output][idx++][0] = sigma*e.evaluate();
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate();
    }
    output++;

    // Average temperature Output
    M_betaFqm[output++][0][0] = 1;

    // Intensity Output
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( exAtM.expression() != "0" )
            M_betaFqm[output++][0][0] = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
    }
}

ThermoElectric::beta_vector_type
ThermoElectric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

ThermoElectric::affine_decomposition_type
ThermoElectric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<ThermoElectric::sparse_matrix_ptrtype> >
ThermoElectric::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<ThermoElectric::element_ptrtype> >
ThermoElectric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

ThermoElectric::element_type
ThermoElectric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto V = Vh->element();
    auto phiV = Vh->element();
    auto T = Th->element();
    auto phiT = Th->element();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    auto gamma = doption("thermoelectric.gamma");

    /***************************** Electro *****************************/
    tic();
    auto aV = form2(_test=Vh, _trial=Vh);
    // V
    for( auto const& mat : materials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        aV += integrate( markedelements(M_mesh, mat.first),
                         sigma*inner(gradt(V),grad(phiV)) );
    }

    // V Dirichlet condition
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
        aV += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         sigma*(gamma/hFace()*inner(idt(V),id(phiV))
                                -inner(gradt(V)*N(),id(phiV))
                                -inner(grad(phiV)*N(),idt(V)) ) );
    }

    auto fV = form1(_test=Vh);
    for( auto const& exAtM : bc["potential"]["SourceTerm"] )
    {
        auto e1 = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e1.expression().hasSymbol(param.first) )
                e1.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto e2 = expr(exAtM.expression2());

        fV += integrate( markedelements(M_mesh, exAtM.marker() ),
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

        fV += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         sigma*e.evaluate()*(gamma/hFace()*id(phiV) -  grad(phiV)*N()) );
    }
    aV.solve( _solution=V, _rhs=fV, _name="electro" );
    toc("electro");

    /***************************** Thermo *****************************/
    tic();
    auto aT = form2(_test=Th, _trial=Th);
    // T
    for( auto const& mat : materials )
    {
        auto k = mu.parameterNamed(mat.second.getString("kKey"));
        aT += integrate( markedelements(M_mesh, mat.first),
                         k*inner(gradt(T),grad(phiT)) );
    }

    // T Robin Condition
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        aT += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         e.evaluate()*inner(idt(T), id(phiT)) );
    }

    auto fT = form1(_test=Th);
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e1 = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e1.expression().hasSymbol(param.first) )
                e1.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto e2 = expr(exAtM.expression2());

        fT += integrate( markedelements(M_mesh, exAtM.marker() ),
                         e2*e1.evaluate()*id(phiT) );
    }
    // T Robin condition
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        fT += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         e.evaluate()*id(phiT) );
    }
    aT.solve( _solution=T, _rhs=fT, _name="thermo" );
    toc("thermo");

    if( boption("thermoelectric.export-FE") )
    {
        auto e = exporter(M_mesh);
        auto j = M_Jh->element();
        for( auto const& mat : materials )
        {
            auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
            j += vf::project( _range=markedelements(M_mesh, mat.first), _space=M_Jh,
                              _expr=sigma*inner(gradv(V),gradv(V)) );
        }
        e->add( "joule", j );
        e->add("V", V);
        e->add("T", T);
        e->save();
    }

    auto solution = Xh->element();
    solution.template element<0>() = V;
    solution.template element<1>() = T;

    return solution;
}

double ThermoElectric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index < Nl() )
    {
        for ( int q = 0; q < Ql(output_index); q++ )
        {
            for( int m = 0; m < mLQF(output_index, q); ++m )
            {
                element_ptrtype eltF( new element_type( Xh ) );
                *eltF = *M_Fqm[output_index][q][0];
                output += M_betaFqm[output_index][q][0]*dot( *eltF, u );
                //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
            }
        }
    }
    else
        throw std::logic_error( "[ThermoElectric::output] error with output_index : only 0, 1, 2 " );
    return output;
}

FEELPP_CRB_PLUGIN( ThermoElectric, thermoelectricsimple)
}
