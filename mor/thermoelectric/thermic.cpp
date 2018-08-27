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
#include "thermic.hpp"

namespace Feel {

Thermic::Thermic()
    : super_type( "thermic" )
{}

Thermic::Thermic( mesh_ptrtype mesh )
    : super_type( "thermic" )
{
    this->M_mesh = mesh;
}

int Thermic::Qa()
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials();
    int robSize = bc["temperature"]["Robin"].size();
    int matSize = materials.materialWithPhysic("thermic").size();
    return matSize+robSize;
}

int Thermic::mQA( int q )
{
    return 1;
}

int Thermic::Nl()
{
    return 2;
}

int Thermic::Ql( int l)
{
    auto bc = M_modelProps->boundaryConditions();
    int sourceSize = bc["temperature"]["SourceTerm"].size();
    int robSize = bc["temperature"]["Robin"].size();
    switch(l) {
    case 0:
        return sourceSize+robSize;
    case 1:
        return 1;
    default:
        return 0;
    }
}

int Thermic::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    case 1:
        return mAverageTempQ(q);
    default:
        return 0;
    }
}

int Thermic::mCompliantQ(int q )
{
    return 1;
}

int Thermic::mAverageTempQ(int q )
{
    return 1;
}

void Thermic::resizeQm( bool resizeMatrix )
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

Thermic::parameter_type
Thermic::paramFromProperties() const
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

void Thermic::initModel()
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
    auto thermoMat = materials.materialWithPhysic("thermic");
    for( auto const& mat : thermoMat )
        range.push_back(mat.first);

    auto domain = markedelements(M_mesh, range);
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=domain ) );

    Feel::cout << "Temperature nDof  : " << Xh->nDof() << std::endl;

    if( !M_T )
        M_T = element_ptrtype( new element_type( Xh ) );

    this->resizeQm();
    this->decomposition();
}

void Thermic::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
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

void Thermic::decomposition()
{
    auto T = Xh->element();
    auto phiT = Xh->element();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermic");

    /************** Left hand side **************/
    // thermo
    int idx = 0;
    for( auto const& mat : materials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(T),grad(phiT)) );
        M_Aqm[idx++][0] = a0.matrixPtr();
    }

    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto aTR = form2(_test=Xh, _trial=Xh);
        aTR += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          inner(idt(T),id(phiT)) );
        M_Aqm[idx++][0] = aTR.matrixPtr();
    }

    /************** Right hand side **************/
    idx = 0;
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression2());
        auto fTST = form1(_test=Xh);
        fTST = integrate( markedelements(M_mesh, exAtM.marker() ),
                          e*id(phiT) );
        M_Fqm[0][idx++][0] = fTST.vectorPtr();
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto fTR = form1(_test=Xh);
        fTR = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         id(phiT) );
        M_Fqm[0][idx++][0] = fTR.vectorPtr();
    }

    // Average temperature output
    auto fAvgT = form1(_test=Xh);
    for( auto const& mat : materials )
    {
        double area = integrate( markedelements(M_mesh, mat.first), cst(1.0) ).evaluate()(0,0) ;
        fAvgT += integrate( markedelements(M_mesh, mat.first ),
                            id(phiT)/cst(area) );
    }
    M_Fqm[1][0][0] = fAvgT.vectorPtr();

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(gradt(T),grad(phiT)) );
    M_energy_matrix = m.matrixPtr();
}

Thermic::beta_vector_type
Thermic::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
}

Thermic::beta_type
Thermic::computeBetaQm( element_type const& T, parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

Thermic::beta_type
Thermic::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

Thermic::beta_type
Thermic::computeBetaQm( parameter_type const& mu )
{
    this->fillBetaQm(mu);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void Thermic::fillBetaQm( parameter_type const& mu )
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermic");

    // Left Hand Side
    int idx = 0;
    for( auto const& mat : materials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("kKey"));
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaAqm[idx++][0] = e.evaluate();
    }

    // Right Hand Side
    idx = 0;
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[0][idx++][0] = e.evaluate();
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[0][idx++][0] = e.evaluate();
    }

    // Average temperature Output
    M_betaFqm[1][0][0] = 1;
}

Thermic::beta_vector_type
Thermic::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

Thermic::affine_decomposition_type
Thermic::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<Thermic::sparse_matrix_ptrtype> >
Thermic::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<Thermic::element_ptrtype> >
Thermic::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

Thermic::element_type
Thermic::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto T = Xh->element();
    auto phiT = Xh->element();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermic");

    /***************************** Thermo *****************************/
    tic();
    auto a = form2(_test=Xh, _trial=Xh);
    // T
    for( auto const& mat : materials )
    {
        auto k = mu.parameterNamed(mat.second.getString("kKey"));
        a += integrate( markedelements(M_mesh, mat.first),
                        k*inner(gradt(T),grad(phiT)) );
    }

    // T Robin Condition
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                        e.evaluate()*inner(idt(T), id(phiT)) );
    }

    auto f = form1(_test=Xh);
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e1 = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e1.expression().hasSymbol(param.first) )
                e1.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto e2 = expr(exAtM.expression2());

        f += integrate( markedelements(M_mesh, exAtM.marker() ),
                        e2*e1.evaluate()*id(phiT) );
    }
    // T Robin condition
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                        e.evaluate()*id(phiT) );
    }
    a.solve( _solution=T, _rhs=f, _name="thermo" );

    if( boption("thermoelectric.export-FE") )
    {
        auto e = exporter(M_mesh);
        e->add("T", T);
        e->save();
    }
    toc("mono");

    return T;
}

double Thermic::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
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

FEELPP_CRB_PLUGIN( Thermic, thermic)
}
