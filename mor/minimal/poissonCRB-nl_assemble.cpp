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

#include "poissonCRB-nl.hpp"

namespace Feel {

void
PoissonNL::assemble()
{
    this->resize();

    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto gamma = doption("poisson.gamma");
    auto bdConditions = M_modelProps->boundaryConditions2();

    /**********  left hand side *********/
    int i = 0;
    for( auto const& [key,mat] : M_therMaterials )
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


    for( auto const& [key,mat] : M_elecMaterials )
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
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        int index = this->indexOfElecMat(bd.material());
        auto eim_sigma = this->scalarContinuousEim()[M_nbTherMat+index];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto a2 = form2( _test=Xh, _trial=Xh );
            a2 += integrate( markedfaces(M_mesh,bd.markers()),
                             idv(eim_sigma->q(m))*(gamma/hFace()*inner(idt(u1),id(u1))
                                                   -inner(gradt(u1)*N(),id(u1))
                                                   -inner(grad(u1)*N(),idt(u1)) ) );
            M_Aqm[i][m] = a2.matrixPtr();
        }
        ++i;
    }

    // Robin
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");
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
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto eim_sigma = this->scalarContinuousEim()[M_nbTherMat+i];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            for( int mm = 0; mm < eim_grad->mMax(); ++mm )
            {
                auto f0 = form1(_test=Xh);
                f0 = integrate( markedelements(M_mesh, mat.meshMarkers()),
                                idv(eim_grad->q(mm))*idv(eim_sigma->q(m))*id(u2) );
                M_Fqm[0][i][mm+m*eim_grad->mMax()] = f0.vectorPtr();
            }
        }
        ++i;
    }

    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        int index = this->indexOfElecMat(bd.material());
        auto eim_sigma = this->scalarContinuousEim()[M_nbTherMat+index];
        for( int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto f1 = form1(_test=Xh);
            f1 = integrate( markedfaces(M_mesh, bd.markers()),
                            idv(eim_sigma->q(m))*(gamma/hFace()*id(u1) - grad(u1)*N() ) );
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
            auto eimSigma = this->scalarContinuousEim()[indexOfElecMat(mat)+M_nbTherMat];
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

PoissonNL::betaqm_type
PoissonNL::computeBetaQm( parameter_type const& mu )
{
    std::vector<vectorN_type> betaEimK(M_nbTherMat);
    for( int i = 0; i < M_nbTherMat; ++i )
        betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu);
    std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
    for( int i = 0; i < M_nbElecMat; ++i )
        betaEimSigma[i] = this->scalarContinuousEim()[M_nbTherMat+i]->beta(mu);
    auto betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu);
    this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

PoissonNL::betaqm_type
PoissonNL::computeBetaQm( element_type const& u, parameter_type const& mu)
{
    std::vector<vectorN_type> betaEimK(M_nbTherMat);
    for( int i = 0; i < M_nbTherMat; ++i )
        betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu, u);
    std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
    for( int i = 0; i < M_nbElecMat; ++i )
        betaEimSigma[i] = this->scalarContinuousEim()[M_nbTherMat+i]->beta(mu, u);
    auto betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu, u);
    this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

PoissonNL::betaqm_type
PoissonNL::computeBetaQm( vectorN_type const& urb, parameter_type const& mu)
{
    std::vector<vectorN_type> betaEimK(M_nbTherMat);
    for( int i = 0; i < M_nbTherMat; ++i )
        betaEimK[i] = this->scalarContinuousEim()[i]->beta(mu, urb);
    std::vector<vectorN_type> betaEimSigma(M_nbElecMat);
    for( int i = 0; i < M_nbElecMat; ++i )
        betaEimSigma[i] = this->scalarContinuousEim()[M_nbTherMat+i]->beta(mu, urb);
    auto betaEimGrad = this->scalarDiscontinuousEim()[0]->beta(mu, urb);
    this->fillBetaQm(mu, betaEimGrad, betaEimK, betaEimSigma);

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

void PoissonNL::fillBetaQm( parameter_type const& mu, vectorN_type const& betaGrad, std::vector<vectorN_type> const& betaK, std::vector<vectorN_type> const& betaSigma )
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    /**********  left hand side *********/
    int i = 0;
    for( auto const& mat : M_therMaterials )
    {
        for( int m = 0; m < betaK[i].size(); ++m )
            M_betaAqm[i][m] = betaK[i](m);
        ++i;
    }
    for( auto const& mat : M_elecMaterials )
    {
        for( int m = 0; m < betaSigma[i-M_nbTherMat].size(); ++m )
            M_betaAqm[i][m] = betaSigma[i-M_nbTherMat](m);
        ++i;
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        for( int m = 0; m < betaSigma[0].size(); ++m )
            M_betaAqm[i][m] = betaSigma[0](m);
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
    for( auto const& mat : M_elecMaterials )
    {
        for( int m = 0; m < betaSigma[i].size(); ++m )
            for( int mm = 0; mm < betaGrad.size(); ++mm )
                M_betaFqm[0][i][mm+m*betaGrad.size()] = betaGrad(mm)*betaSigma[i](m);
        ++i;
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto e = bd.expr();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        for( int m = 0; m < betaSigma[0].size(); ++m )
            M_betaFqm[0][i][m] = e.evaluate()(0,0)*betaSigma[0](m);
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
        }
        else if( output.type() == "intensity" )
        {
            int index = indexOfElecMat(output.getString("material"));
            for( int m = 0; m < betaSigma[index].size(); ++m )
                M_betaFqm[i][0][m] = betaSigma[index](m);
        }
        ++i;
    }
}


std::vector<std::vector<PoissonNL::element_ptrtype> >
PoissonNL::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

PoissonNL::beta_vector_type
PoissonNL::computeBetaInitialGuess( parameter_type const& mu)
{
    auto mu0 = this->Dmu->min();
    M_betaInitialGuess[0][0] = 1;
    M_betaInitialGuess[1][0] = 1;
    return this->M_betaInitialGuess;
}

std::vector< std::vector<PoissonNL::sparse_matrix_ptrtype> >
PoissonNL::computeLinearDecompositionA()
{
    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto gamma = doption("poisson.gamma");
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto aqm = std::vector<std::vector<sparse_matrix_ptrtype> >(Qa(), std::vector<sparse_matrix_ptrtype>(1));

    int i = 0;
    for( auto const& mat : M_therMaterials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.second.meshMarkers()),
                        inner(gradt(u2),grad(u2)) );
        aqm[i][0] = a0.matrixPtr();
        ++i;
    }
    for( auto const& mat : M_elecMaterials )
    {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 = integrate( markedelements(M_mesh, mat.second.meshMarkers()),
                        inner(gradt(u1),grad(u1)) );
        aqm[i][0] = a1.matrixPtr();
        ++i;
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto a2 = form2(_test=Xh, _trial=Xh);
        a2 += integrate( markedfaces(M_mesh, bd.markers()),
                         gamma/hFace()*inner(idt(u1),id(u1))
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
        i++;
    }
    return aqm;
}

PoissonNL::beta_vector_type
PoissonNL::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    beta_vector_type beta(Qa(), std::vector<double>(1));

    int i = 0;
    for( auto const& [key,mat] : M_therMaterials )
    {
        beta[i][0] = mu.parameterNamed(mat.getString("misc.sigmaKey"))*mu.parameterNamed("L")*293;
        ++i;
    }
    for( auto const& [key,mat] : M_elecMaterials )
    {
        beta[i][0] = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        ++i;
    }
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        beta[i][0] = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        ++i;
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        beta[i][0] = e.evaluate()(0,0);
        ++i;
    }
    return beta;
}

PoissonNL::element_type
PoissonNL::solve(parameter_type const& mu)
{
    tic();
    auto gamma = doption("poisson.gamma");

    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto uOld1 = M_InitialGuess[0][0]->template element<0>();
    auto uOld2 = M_InitialGuess[1][0]->template element<1>();

    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    double tol = doption("poisson.tolerance");
    int maxit = ioption("poisson.maxit");
    double increment = 0;
    int it = 0;

    do
    {
        auto a = form2(_test=Xh, _trial=Xh);
        for( auto const& [key,mat] : M_therMaterials )
        {
            auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
            auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
            auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
            auto k = sigma*mu.parameterNamed("L")*idv(uOld2);
            a += integrate( markedelements(M_mesh, mat.meshMarkers()),
                            k*inner(gradt(u2),grad(u2)) );
        }
        for( auto const& [key,mat] : M_elecMaterials )
        {
            auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
            auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
            auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
            a += integrate( markedelements(M_mesh, mat.meshMarkers()),
                            sigma*inner(gradt(u1),grad(u1)) );
        }

        // Dirichlet condition
        for( auto const&[bdName, bd] : potentialDirichlet )
        {
            auto mat = M_elecMaterials[bd.material()];
            auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
            auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
            auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
            a += integrate( markedfaces(M_mesh,bd.markers()),
                            sigma*(gamma/hFace()*inner(idt(u1),id(u1))
                                   -inner(gradt(u1)*N(),id(u1))
                                   -inner(grad(u1)*N(),idt(u1)) ) );
        }

        for( auto const&[bdName, bd] : temperatureRobin )
        {
            auto e = bd.expr1();
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
            a += integrate( markedfaces(M_mesh,bd.markers()),
                            e.evaluate()(0,0)*inner(idt(u2),id(u2)) );
        }

        auto f = form1(_test=Xh);
        for( auto const& [key,mat] : M_elecMaterials )
        {
            auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
            auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
            auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
            f += integrate( markedelements(M_mesh,mat.meshMarkers()),
                            sigma*inner(gradv(uOld1))*id(u2) );
        }
        // Dirichlet condition
        for( auto const&[bdName, bd] : potentialDirichlet )
        {
            auto mat = M_elecMaterials[bd.material()];
            auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
            auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
            auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));

            auto e = bd.expr();
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

            f += integrate( markedfaces(M_mesh, bd.markers()),
                           e.evaluate()(0,0)*sigma*(gamma/hFace()*id(u1) - grad(u1)*N() ) );
        }
        for( auto const&[bdName, bd] : temperatureRobin )
        {
            auto e = bd.expr2();
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
            f += integrate( markedfaces(M_mesh, bd.markers()),
                            e.evaluate()(0,0)*id(u2) );
        }

        a.solve(_solution=U, _rhs=f, _name="thermo-electro");

        increment = normL2(_range=elements(M_mesh), _expr=idv(u2)-idv(uOld2) + idv(u1) - idv(uOld1));
        uOld1 = u1;
        uOld2 = u2;
        ++it;
    }
    while( increment > tol && it < maxit );

    toc("solve");

    return U;
}

PoissonNL::element_type
PoissonNL::solveLinear(parameter_type const& mu)
{
    tic();
    auto gamma = doption("poisson.gamma");
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto XhV = Xh->template functionSpace<0>();
    auto XhT = Xh->template functionSpace<1>();
    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();

    auto aV = form2(_test=XhV, _trial=XhV);
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        aV += integrate( markedelements(M_mesh, mat.meshMarkers()),
                        sigma0*inner(gradt(u1),grad(u1)) );
    }
    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        aV += integrate( markedfaces(M_mesh,bd.markers()),
                         sigma0*(gamma/hFace()*inner(idt(u1),id(u1))
                                 -inner(gradt(u1)*N(),id(u1))
                                 -inner(grad(u1)*N(),idt(u1)) ) );
    }

    // Dirichlet condition
    auto fV = form1(_test=XhV);
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));

        auto e = bd.expr();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        fV += integrate( markedfaces(M_mesh, bd.markers()),
                        e.evaluate()(0,0)*sigma0*(gamma/hFace()*id(u1) - grad(u1)*N() ) );
    }

    aV.solve(_solution=u1, _rhs=fV, _name="electro");

    auto aT = form2(_test=XhT, _trial=XhT);
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        auto k = sigma0*mu.parameterNamed("L")*293;
        aT += integrate( markedelements(M_mesh,mat.meshMarkers()), k*inner(gradt(u2),grad(u2)) );
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        aT += integrate( markedfaces(M_mesh,bd.markers()),
                         e.evaluate()(0,0)*inner(idt(u2),id(u2)) );
    }

    auto fT = form1(_test=XhT);
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        fT += integrate( markedelements(M_mesh, mat.meshMarkers()),
                         sigma0*inner(gradv(u1))*id(u2) );
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        fT += integrate( markedfaces(M_mesh, bd.markers()),
                         e.evaluate()(0,0)*id(u2) );
    }

    aT.solve(_solution=u2, _rhs=fT, _name="thermo");

    toc("solveLinear");

    return U;
}

}
