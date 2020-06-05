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

#include "thermoelectric-nl.hpp"

namespace Feel {

void
ThermoElectricNL::assemble()
{
    this->resize();

    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    int i = 0;
    /**********  left hand side *********/
    for( auto const& [name,mat] : M_therMaterials )
    {
        auto eim_k = this->scalarContinuousEim()[i];
        for( int m = 0; m < eim_k->mMax(); ++m )
        {
            auto a0 = form2(_test=Xh, _trial=Xh);
            a0 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                            idv(eim_k->q(m))*inner(gradt(u2),grad(u2)) );
            M_Aqm[i][m] = a0.matrixPtr();
        }
        ++i;
    }


    for( auto const& [name,mat] : M_elecMaterials )
    {
        auto eim_sigma = this->scalarContinuousEim()[i];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto a1 = form2( _test=Xh, _trial=Xh );
            a1 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                            idv(eim_sigma->q(m))*inner(gradt(u1),grad(u1)) );
            M_Aqm[i][m] = a1.matrixPtr();
        }
        ++i;
    }

    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[bd.material()]];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto a2 = form2( _test=Xh, _trial=Xh );
            a2 = integrate( markedfaces(M_mesh,bd.markers()),
                            idv(eim_sigma->q(m))*(M_gamma/hFace()*inner(idt(u1),id(u1))
                                                  -inner(gradt(u1)*N(),id(u1))
                                                  -inner(grad(u1)*N(),idt(u1)) ) );
            M_Aqm[i][m] = a2.matrixPtr();
        }
        ++i;
    }

    // Robin
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto a3 = form2(_test=Xh, _trial=Xh);
        a3 = integrate( markedfaces(M_mesh, bd.markers()), inner(idt(u2),id(u2)) );
        M_Aqm[i][0] = a3.matrixPtr();
        ++i;
    }

    /**********  right hand side *********/
    i = 0;
    auto eim_grad = this->scalarDiscontinuousEim()[0];
    for( auto const& [name,mat] : M_elecMaterials )
    {
        auto eim_sigma = this->scalarContinuousEim()[i+M_nbTherMat];
        int M = 0;
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            for( int mm = 0; mm < eim_grad->mMax(); ++mm )
            {
                auto f0 = form1(_test=Xh);
                f0 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                                M_sigmaMax*idv(eim_sigma->q(m))*idv(eim_grad->q(mm))*id(u2) );
                M_Fqm[0][i][M++] = f0.vectorPtr();
            }
        }
        ++i;
    }

    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[bd.material()]];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto f1 = form1(_test=Xh);
            f1 = integrate( markedfaces(M_mesh, bd.markers()),
                            idv(eim_sigma->q(m))*(M_gamma/hFace()*id(u1) - grad(u1)*N() ) );
            M_Fqm[0][i][m] = f1.vectorPtr();
        }
        ++i;
    }

    // Robin
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto f2 = form1(_test=Xh);
        f2 = integrate( markedfaces(M_mesh,bd.markers()), id(u2) );
        M_Fqm[0][i][0] = f2.vectorPtr();
        ++i;
    }

    /********* Outputs ********/
    i = 1;
    auto outputs = M_modelProps->outputs();
    for( auto const& [key, output] : outputs )
    {
        if( output.type() == "averageTemp" )
        {
            auto dim = output.dim();
            auto fAvgT = form1(_test=Xh);
            if( dim == 3 )
            {
                auto range = markedelements(M_mesh, output.markers());
                double area = integrate(_range=range, _expr=cst(1.)).evaluate()(0,0);
                fAvgT = integrate( _range=range, _expr=id(u2)/cst(area) );
            }
            else if( dim == 2 )
            {
                auto range = markedfaces(M_mesh, output.markers() );
                double area = integrate(_range=range, _expr=cst(1.)).evaluate()(0,0);
                fAvgT = integrate( _range=range, _expr=id(u2)/cst(area) );
            }
            M_Fqm[i++][0][0] = fAvgT.vectorPtr();
        }
        else if( output.type() == "intensity" )
        {
            auto mat = output.getString("material");
            auto eimSigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
            for( int m = 0; m < eimSigma->mMax(); ++m )
            {
                auto fI = form1(_test=Xh);
                fI = integrate( markedfaces(M_mesh,output.markers() ),
                                -idv(eimSigma->q(m))*grad(u1)*N() );
                M_Fqm[i][0][m] = fI.vectorPtr();
            }
            ++i;
        }
    }
}

ThermoElectricNL::betaqm_type
ThermoElectricNL::computeBetaQm( parameter_type const& mu )
{
    std::vector<vectorN_type> betaEimK(M_nbTherMat);
    for( int i = 0; i < M_nbTherMat; ++i )
        betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu);
    std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
    for( int i = 0; i < M_nbElecMat; ++i )
        betaEimSigma[i] = this->scalarContinuousEim()[i+M_nbTherMat]->beta(mu);
    vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu);
    this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

ThermoElectricNL::betaqm_type
ThermoElectricNL::computeBetaQm( element_type const& u, parameter_type const& mu)
{
    int i = 0;
    std::vector<vectorN_type> betaEimK(M_nbTherMat);
    for( int i = 0; i < M_nbTherMat; ++i )
        betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu, u);
    std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
    for( int i = 0; i < M_nbElecMat; ++i )
        betaEimSigma[i] = this->scalarContinuousEim()[i+M_nbTherMat]->beta(mu, u);
    vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu, u);
    this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

ThermoElectricNL::betaqm_type
ThermoElectricNL::computeBetaQm( vectorN_type const& urb, parameter_type const& mu)
{
    int i = 0;
    std::vector<vectorN_type> betaEimK(M_nbTherMat);
    for( int i = 0; i < M_nbTherMat; ++i )
        betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu, urb);
    std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
    for( int i = 0; i < M_nbElecMat; ++i )
        betaEimSigma[i] = this->scalarContinuousEim()[i+M_nbTherMat]->beta(mu, urb);
    vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu, urb);
    this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

void ThermoElectricNL::fillBetaQm( parameter_type const& mu, vectorN_type const& betaGrad, std::vector<vectorN_type> const& betaK, std::vector<vectorN_type> const& betaSigma )
{
    // Feel::cout << "betaGrad:\n" << betaGrad << std::endl;
    // Feel::cout << "betaK:\n" << betaK << std::endl;
    // Feel::cout << "betaSigma:\n" << betaSigma << std::endl;
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");
    /**********  left hand side *********/
    int i = 0;
    for( i = 0; i < M_nbTherMat; ++i )
        for( int m = 0; m < M_betaAqm[i].size(); ++m )
            M_betaAqm[i][m] = betaK[i](m);
    for( i = M_nbTherMat; i < M_nbTherMat+M_nbElecMat; ++i )
        for( int m = 0; m < M_betaAqm[i].size(); ++m )
            M_betaAqm[i][m] = betaSigma[i-M_nbTherMat](m);
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        int idx = M_elecEimIndex[bd.material()]-M_nbTherMat;
        for( int m = 0; m < M_betaAqm[i].size(); ++m )
            M_betaAqm[i][m] = betaSigma[idx](m);
        ++i;
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        M_betaAqm[i][0] = e.evaluate()(0,0);
        ++i;
    }

    /**********  right hand side *********/
    i = 0;
    for( i = 0; i < M_nbElecMat; ++i )
    {
        int M = 0;
        for( int m = 0; m < M_betaFqm[0][i].size()/M_eimGradSize; ++m )
            for( int mm = 0; mm < M_eimGradSize; ++mm )
                M_betaFqm[0][i][M++] = betaSigma[i](m)*betaGrad(mm);
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto e = bd.expr();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        int idx = M_elecEimIndex[bd.material()]-M_nbTherMat;
        for( int m = 0; m < M_betaFqm[0][i].size(); ++m )
            M_betaFqm[0][i][m] = e.evaluate()(0,0)*betaSigma[idx](m);
        ++i;
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr2();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        M_betaFqm[0][i][0] = e.evaluate()(0,0);
        ++i;
    }

    /**********  outputs *********/
    i = 1;
    auto outputs = M_modelProps->outputs();
    for( auto const& [key, output] : outputs )
    {
        if( output.type() == "averageTemp" )
        {
            M_betaFqm[i][0][0] = 1;
            ++i;
        }
        else if( output.type() == "intensity" )
        {
            auto mat = output.getString("material");
            int idx = M_elecEimIndex[mat] - M_nbTherMat;
            for( int m = 0; m < M_betaFqm[i][0].size(); ++m )
                M_betaFqm[i][0][m] = betaSigma[idx](m);
            ++i;
        }
    }
}


std::vector<std::vector<ThermoElectricNL::element_ptrtype> >
ThermoElectricNL::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

ThermoElectricNL::beta_vector_type
ThermoElectricNL::computeBetaInitialGuess( parameter_type const& mu)
{
    M_betaInitialGuess[0][0] = 1;
    M_betaInitialGuess[1][0] = 1;
    return this->M_betaInitialGuess;
}

std::vector< std::vector<ThermoElectricNL::sparse_matrix_ptrtype> >
ThermoElectricNL::computeLinearDecompositionA()
{
    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();

    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto aqm = std::vector<std::vector<sparse_matrix_ptrtype> >(Qa(), std::vector<sparse_matrix_ptrtype>(1));

    int i = 0;
    for( auto const& [name,mat] : M_therMaterials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh,mat.meshMarkers()),
                        inner(gradt(u2),grad(u2)) );
        aqm[i][0] = a0.matrixPtr();
        ++i;
    }
    for( auto const& [name,mat] : M_elecMaterials )
    {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                        inner(gradt(u1),grad(u1)) );
        aqm[i][0] = a1.matrixPtr();
        ++i;
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto a2 = form2(_test=Xh, _trial=Xh);
        a2 = integrate( markedfaces(M_mesh, bd.markers()),
                        M_gamma/hFace()*inner(idt(u1),id(u1))
                        -inner(gradt(u1)*N(),id(u1))
                        -inner(grad(u1)*N(),idt(u1)) );
        aqm[i][0] = a2.matrixPtr();
        ++i;
    }

    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto a3 = form2(_test=Xh, _trial=Xh);
        a3 = integrate( markedfaces(M_mesh, bd.markers()), inner(idt(u2),id(u2)) );
        aqm[i][0] = a3.matrixPtr();
        ++i;
    }

    return aqm;
}

ThermoElectricNL::beta_vector_type
ThermoElectricNL::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta(Qa(), std::vector<double>(1));
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    int i = 0;
    for( auto const& [name,mat] : M_elecMaterials )
        beta[i++][0] = mu.parameterNamed(mat.getString("misc.sigmaKey"))*mu.parameterNamed("L")*293;
    for( auto const& [name,mat] : M_elecMaterials )
        beta[i++][0] = mu.parameterNamed(mat.getString("misc.sigmaKey"));
    for( auto const&[bdName, bd] : potentialDirichlet )
        beta[i++][0] = mu.parameterNamed(M_elecMaterials[bd.material()].getString("misc.sigmaKey"));
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        beta[i++][0] = e.evaluate()(0,0);
    }
    return beta;
}

}

