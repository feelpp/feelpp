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

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::assemble()
{
    if ( M_useDEIM )
        this->assembleWithDEIM();
    else
        this->assembleWithEIM();
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::assembleWithEIM()
{
    this->resize();

    auto U = this->Xh->element();
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
            auto a0 = form2(_test=this->Xh, _trial=this->Xh);
            a0 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                            idv(eim_k->q(m))*inner(gradt(u2),grad(u2)) );
            this->M_Aqm[i][m] = a0.matrixPtr();
        }
        ++i;
    }


    for( auto const& [name,mat] : M_elecMaterials )
    {
        auto eim_sigma = this->scalarContinuousEim()[i];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto a1 = form2( _test=this->Xh, _trial=this->Xh );
            a1 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                            idv(eim_sigma->q(m))*inner(gradt(u1),grad(u1)) );
            this->M_Aqm[i][m] = a1.matrixPtr();
        }
        ++i;
    }

    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[bd.material()]];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto a2 = form2( _test=this->Xh, _trial=this->Xh );
            a2 = integrate( markedfaces(M_mesh,bd.markers()),
                            idv(eim_sigma->q(m))*(M_gamma/hFace()*inner(idt(u1),id(u1))
                                                  -inner(gradt(u1)*N(),id(u1))
                                                  -inner(grad(u1)*N(),idt(u1)) ) );
            this->M_Aqm[i][m] = a2.matrixPtr();
        }
        ++i;
    }

    // Robin
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto a3 = form2(_test=this->Xh, _trial=this->Xh);
        a3 = integrate( markedfaces(M_mesh, bd.markers()), inner(idt(u2),id(u2)) );
        this->M_Aqm[i][0] = a3.matrixPtr();
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
                auto f0 = form1(_test=this->Xh);
                f0 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                                M_sigmaMax*idv(eim_sigma->q(m))*idv(eim_grad->q(mm))*id(u2) );
                this->M_Fqm[0][i][M++] = f0.vectorPtr();
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
            auto f1 = form1(_test=this->Xh);
            f1 = integrate( markedfaces(M_mesh, bd.markers()),
                            idv(eim_sigma->q(m))*(M_gamma/hFace()*id(u1) - grad(u1)*N() ) );
            this->M_Fqm[0][i][m] = f1.vectorPtr();
        }
        ++i;
    }

    // Robin
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto f2 = form1(_test=this->Xh);
        f2 = integrate( markedfaces(M_mesh,bd.markers()), id(u2) );
        this->M_Fqm[0][i][0] = f2.vectorPtr();
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
            auto fAvgT = form1(_test=this->Xh);
            if( dim == Dim )
            {
                auto range = markedelements(M_mesh, output.markers());
                double area = integrate(_range=range, _expr=cst(1.)).evaluate()(0,0);
                fAvgT = integrate( _range=range, _expr=id(u2)/cst(area) );
            }
            else if( dim == Dim-1 )
            {
                auto range = markedfaces(M_mesh, output.markers() );
                double area = integrate(_range=range, _expr=cst(1.)).evaluate()(0,0);
                fAvgT = integrate( _range=range, _expr=id(u2)/cst(area) );
            }
            this->M_Fqm[i++][0][0] = fAvgT.vectorPtr();
        }
        else if( output.type() == "intensity" )
        {
            auto mat = output.getString("material");
            auto eimSigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
            for( int m = 0; m < eimSigma->mMax(); ++m )
            {
                auto fI = form1(_test=this->Xh);
                fI = integrate( markedfaces(M_mesh,output.markers() ),
                                -idv(eimSigma->q(m))*grad(u1)*N() );
                this->M_Fqm[i][0][m] = fI.vectorPtr();
            }
            ++i;
        }
    }
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::betaqm_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeBetaQm( parameter_type const& mu )
{
    if( M_useDEIM )
    {
        auto betaDEIM = this->deim()->beta(mu);
        auto betaMDEIM = this->mdeim()->beta(mu);
        this->fillBetaQmWithDEIM(mu, betaDEIM, betaMDEIM);
    }
    else
    {
        std::vector<vectorN_type> betaEimK(M_nbTherMat);
        for( int i = 0; i < M_nbTherMat; ++i )
            betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu);
        std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
        for( int i = 0; i < M_nbElecMat; ++i )
            betaEimSigma[i] = this->scalarContinuousEim()[i+M_nbTherMat]->beta(mu);
        vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu);
        this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);
    }

    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm );
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::betaqm_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeBetaQm( element_type const& u, parameter_type const& mu)
{
    if( M_useDEIM )
    {
        auto betaDEIM = this->deim()->beta(mu, u);
        auto betaMDEIM = this->mdeim()->beta(mu, u);
        this->fillBetaQmWithDEIM(mu, betaDEIM, betaMDEIM);
    }
    else
    {
        std::vector<vectorN_type> betaEimK(M_nbTherMat);
        for( int i = 0; i < M_nbTherMat; ++i )
            betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu, u);
        std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
        for( int i = 0; i < M_nbElecMat; ++i )
            betaEimSigma[i] = this->scalarContinuousEim()[i+M_nbTherMat]->beta(mu, u);
        Feel::cout << "using betaEimGrad(u, mu)" << std::endl;
        // vectorN_type betaEimGrad (M_eimGradSize);
        vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu, u);
        this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);
    }

    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm );
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::betaqm_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeBetaQm( vectorN_type const& urb, parameter_type const& mu)
{
    if( M_useDEIM )
    {
        auto betaDEIM = this->deim()->beta(mu, urb);
        auto betaMDEIM = this->mdeim()->beta(mu, urb);
        this->fillBetaQmWithDEIM(mu, betaDEIM, betaMDEIM);
    }
    else
    {
        std::vector<vectorN_type> betaEimK(M_nbTherMat);
        for( int i = 0; i < M_nbTherMat; ++i )
            betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu, urb);
        std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
        for( int i = 0; i < M_nbElecMat; ++i )
            betaEimSigma[i] = this->scalarContinuousEim()[i+M_nbTherMat]->beta(mu, urb);
        vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu, urb);
        this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);
    }

    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm );
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::fillBetaQm( parameter_type const& mu,
                                                   vectorN_type const& betaGrad,
                                                   std::vector<vectorN_type> const& betaK,
                                                   std::vector<vectorN_type> const& betaSigma )
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
        for( int m = 0; m < this->M_betaAqm[i].size(); ++m )
            this->M_betaAqm[i][m] = betaK[i](m);
    for( i = M_nbTherMat; i < M_nbTherMat+M_nbElecMat; ++i )
        for( int m = 0; m < this->M_betaAqm[i].size(); ++m )
            this->M_betaAqm[i][m] = betaSigma[i-M_nbTherMat](m);
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        int idx = M_elecEimIndex[bd.material()]-M_nbTherMat;
        for( int m = 0; m < this->M_betaAqm[i].size(); ++m )
            this->M_betaAqm[i][m] = betaSigma[idx](m);
        ++i;
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        this->M_betaAqm[i][0] = e.evaluate()(0,0);
        ++i;
    }

    /**********  right hand side *********/
    i = 0;
    for( i = 0; i < M_nbElecMat; ++i )
    {
        int M = 0;
        for( int m = 0; m < this->M_betaFqm[0][i].size()/M_eimGradSize; ++m )
            for( int mm = 0; mm < M_eimGradSize; ++mm )
                this->M_betaFqm[0][i][M++] = betaSigma[i](m)*betaGrad(mm);
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto e = bd.expr();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        int idx = M_elecEimIndex[bd.material()]-M_nbTherMat;
        for( int m = 0; m < this->M_betaFqm[0][i].size(); ++m )
            this->M_betaFqm[0][i][m] = e.evaluate()(0,0)*betaSigma[idx](m);
        ++i;
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr2();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        this->M_betaFqm[0][i][0] = e.evaluate()(0,0);
        ++i;
    }

    /**********  outputs *********/
    i = 1;
    auto outputs = M_modelProps->outputs();
    for( auto const& [key, output] : outputs )
    {
        if( output.type() == "averageTemp" )
        {
            this->M_betaFqm[i][0][0] = 1;
            ++i;
        }
        else if( output.type() == "intensity" )
        {
            auto mat = output.getString("material");
            int idx = M_elecEimIndex[mat] - M_nbTherMat;
            for( int m = 0; m < this->M_betaFqm[i][0].size(); ++m )
                this->M_betaFqm[i][0][m] = betaSigma[idx](m);
            ++i;
        }
    }
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::fillBetaQmWithDEIM( parameter_type const& mu,
                                                           vectorN_type const& betaDEIM,
                                                           vectorN_type const& betaMDEIM )
{
    int A = this->mdeim()->size();
    for( int i = 0; i < A; ++i )
        this->M_betaAqm[0][i] = betaMDEIM(i);

    int F = this->deim()->size();
    for( int i = 0; i < F; ++i )
        this->M_betaFqm[0][0][i] = betaDEIM(i);

}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
std::vector<std::vector<typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::element_ptrtype> >
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::beta_vector_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeBetaInitialGuess( parameter_type const& mu)
{
    this->M_betaInitialGuess[0][0] = 1;
    this->M_betaInitialGuess[1][0] = 1;
    return this->M_betaInitialGuess;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
std::vector< std::vector<typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::sparse_matrix_ptrtype> >
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeLinearDecompositionA()
{
    auto U = this->Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();

    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto aqm = std::vector<std::vector<sparse_matrix_ptrtype> >(Qa(true), std::vector<sparse_matrix_ptrtype>(1));

    int i = 0;
    for( auto const& [name,mat] : M_therMaterials )
    {
        auto a0 = form2(_test=this->Xh, _trial=this->Xh);
        a0 = integrate( markedelements(M_mesh,mat.meshMarkers()),
                        inner(gradt(u2),grad(u2)) );
        aqm[i][0] = a0.matrixPtr();
        ++i;
    }
    for( auto const& [name,mat] : M_elecMaterials )
    {
        auto a1 = form2(_test=this->Xh, _trial=this->Xh);
        a1 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                        inner(gradt(u1),grad(u1)) );
        aqm[i][0] = a1.matrixPtr();
        ++i;
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto a2 = form2(_test=this->Xh, _trial=this->Xh);
        a2 = integrate( markedfaces(M_mesh, bd.markers()),
                        M_gamma/hFace()*inner(idt(u1),id(u1))
                        -inner(gradt(u1)*N(),id(u1))
                        -inner(grad(u1)*N(),idt(u1)) );
        aqm[i][0] = a2.matrixPtr();
        ++i;
    }

    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto a3 = form2(_test=this->Xh, _trial=this->Xh);
        a3 = integrate( markedfaces(M_mesh, bd.markers()), inner(idt(u2),id(u2)) );
        aqm[i][0] = a3.matrixPtr();
        ++i;
    }

    return aqm;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::beta_vector_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta(Qa(true), std::vector<double>(1));
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

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
void
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::assembleWithDEIM()
{
    this->resize();

    int A = this->mdeim()->size();
    auto qa = this->mdeim()->q();
    for( int i = 0; i < A; ++i )
        this->M_Aqm[0][i] = qa[i];

    int F = this->deim()->size();
    auto qf = this->deim()->q();
    for( int i = 0; i < F; ++i )
        this->M_Fqm[0][0][i] = qf[i];

#if 0
    /********* Outputs ********/
    int i = 1;
    auto outputs = M_modelProps->outputs();
    for( auto const& [key, output] : outputs )
    {
        if( output.type() == "averageTemp" )
        {
            auto dim = output.dim();
            auto fAvgT = form1(_test=this->Xh);
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
            this->M_Fqm[i++][0][0] = fAvgT.vectorPtr();
        }
        else if( output.type() == "intensity" )
        {
            auto mat = output.getString("material");
            auto eimSigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
            for( int m = 0; m < eimSigma->mMax(); ++m )
            {
                auto fI = form1(_test=this->Xh);
                fI = integrate( markedfaces(M_mesh,output.markers() ),
                                -idv(eimSigma->q(m))*grad(u1)*N() );
                this->M_Fqm[i][0][m] = fI.vectorPtr();
            }
            ++i;
        }
    }
#endif
}

}

