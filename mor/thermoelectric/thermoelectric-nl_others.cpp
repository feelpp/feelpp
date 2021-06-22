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
#include <iostream>
#include <string>
#include <list>

#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/pchm.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feel.hpp>

namespace Feel {

ThermoElectricNL::sparse_matrix_ptrtype
ThermoElectricNL::assembleForMDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");
    auto parameters = M_modelProps->parameters();

    auto mesh = Xh->mesh();
    auto XhT = Xh->template functionSpace<1>();
    auto u1 = u.template element<0>();
    auto u2 = u.template element<1>();
    auto uOld1 = u.template element<0>();
    auto uOld2 = u.template element<1>();

    auto parameterValues = parameters.toParameterValues();
    for(auto const& [name, param] : parameters )
        if( param.hasMinMax() )
            parameterValues[name] = mu.parameterNamed(name);

    auto a = form2(_test=Xh, _trial=Xh);
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto k1 = mat.getScalar("k", {"heat_T"}, {idv(uOld2)}, parameterValues);
        auto k = XhT->element(k1);
        a += integrate( markedelements(mesh, mat.meshMarkers()),
                        idv(k)*inner(gradt(u2),grad(u2)) );
    }
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma1 = mat.getScalar("sigma", {"heat_T"}, {idv(uOld2)}, parameterValues);
        auto sigma = XhT->element(sigma1);
        a += integrate( markedelements(mesh, mat.meshMarkers()),
                        idv(sigma)*inner(gradt(u1),grad(u1)) );
    }

    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma1 = mat.getScalar("sigma", {"heat_T"}, {idv(uOld2)}, parameterValues);
        auto sigma = XhT->element(sigma1);
        a += integrate( markedfaces(mesh,bd.markers()),
                        idv(sigma)*(M_gamma/hFace()*inner(idt(u1),id(u1))
                                    -inner(gradt(u1)*N(),id(u1))
                                    -inner(grad(u1)*N(),idt(u1)) ) );
    }

    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : parameters )
        {
            if( e.expression().hasSymbol(param.first) )
            {
                if( parameters[param.first].hasMinMax() )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
                else
                    e.setParameterValues( { param.first, parameters[param.first].value() } );
            }
        }
        a += integrate( markedfaces(mesh,bd.markers()),
                        e.evaluate()(0,0)*inner(idt(u2),id(u2)) );
    }

    auto am = a.matrixPtr();
    am->close();

    return am;
}

ThermoElectricNL::vector_ptrtype
ThermoElectricNL::assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");
    auto parameters = M_modelProps->parameters();

    auto mesh = Xh->mesh();
    auto XhT = Xh->template functionSpace<1>();
    auto u1 = u.template element<0>();
    auto u2 = u.template element<1>();
    auto uOld1 = u.template element<0>();
    auto uOld2 = u.template element<1>();

    auto parameterValues = parameters.toParameterValues();
    for(auto const& [name, param] : parameters )
        if( param.hasMinMax() )
            parameterValues[name] = mu.parameterNamed(name);

    auto f = form1(_test=Xh);

    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma1 = mat.getScalar("sigma", {"heat_T"}, {idv(uOld2)}, parameterValues);
        auto sigma = XhT->element(sigma1);
        f += integrate( markedelements(mesh,mat.meshMarkers()),
                        idv(sigma)*inner(gradv(uOld1))*id(u2) );
    }
    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma1 = mat.getScalar("sigma", {"heat_T"}, {idv(uOld2)}, parameterValues);
        auto sigma = XhT->element(sigma1);

        auto e = bd.expr();
        for( auto const& param : parameters )
        {
            if( e.expression().hasSymbol(param.first) )
            {
                if( parameters[param.first].hasMinMax() )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
                else
                    e.setParameterValues( { param.first, parameters[param.first].value() } );
            }
        }

        f += integrate( markedfaces(mesh, bd.markers()),
                        e*idv(sigma)*(M_gamma/hFace()*id(u1) - grad(u1)*N() ) );
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr2();
        for( auto const& param : parameters )
        {
            if( e.expression().hasSymbol(param.first) )
            {
                if( parameters[param.first].hasMinMax() )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
                else
                    e.setParameterValues( { param.first, parameters[param.first].value() } );
            }
        }
        f += integrate( markedfaces(mesh, bd.markers()),
                        e*id(u2) );
    }

    auto fv = f.vectorPtr();
    fv->close();

    return fv;
}

ThermoElectricNL::element_type
ThermoElectricNL::solve(parameter_type const& mu)
{
    tic();

    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto UOld = Xh->element();
    auto uOld1 = UOld.template element<0>();
    auto uOld2 = UOld.template element<1>();
    uOld1 = M_InitialGuess[0][0]->template element<0>();
    uOld2 = M_InitialGuess[1][0]->template element<1>();

    double increment = 0;
    int it = 0;

    do
    {
        auto am = this->assembleForMDEIMnl(mu, UOld, 0);
        auto a = form2(_test=Xh, _trial=Xh, _matrix=am);
        auto fv = this->assembleForDEIMnl(mu, UOld, 0);
        auto f = form1(_test=Xh, _vector=fv);
        a.solve(_solution=U, _rhs=f, _name="thermo-electro");

        increment = normL2(_range=elements(M_mesh), _expr=idv(u2)-idv(uOld2)) / normL2(_range=elements(M_mesh), _expr=idv(uOld2));
        increment += normL2(_range=elements(M_mesh), _expr=idv(u1) - idv(uOld1)) / normL2(_range=elements(M_mesh), _expr=idv(uOld1));
        if( M_verbose > 0 )
            Feel::cout << "iteration " << it << " increment = " << increment << std::endl;

        UOld.template element<0>() = u1;
        UOld.template element<1>() = u2;
        uOld1 = u1;
        uOld2 = u2;
        ++it;
    }
    while( increment > M_tolerance && it < M_maxit );

    toc("solve");

    return U;
}

ThermoElectricNL::element_type
ThermoElectricNL::solveLinear(parameter_type const& mu)
{
    tic();
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto parameters = M_modelProps->parameters();
    auto parameterValues = parameters.toParameterValues();
    for(auto const& [name, param] : parameters )
        if( param.hasMinMax() )
            parameterValues[name] = mu.parameterNamed(name);
    parameterValues["heat_T"] = 293.;

    auto XhV = Xh->template functionSpace<0>();
    auto XhT = Xh->template functionSpace<1>();
    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();

    auto aV = form2(_test=XhV, _trial=XhV);
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma0 = mat.getScalar("sigma", parameterValues);
        aV += integrate( markedelements(M_mesh, mat.meshMarkers()),
                        sigma0*inner(gradt(u1),grad(u1)) );
    }
    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma0 = mat.getScalar("sigma", parameterValues);
        aV += integrate( markedfaces(M_mesh,bd.markers()),
                         sigma0*(M_gamma/hFace()*inner(idt(u1),id(u1))
                                 -inner(gradt(u1)*N(),id(u1))
                                 -inner(grad(u1)*N(),idt(u1)) ) );
    }

    // Dirichlet condition
    auto fV = form1(_test=XhV);
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma0 = mat.getScalar("sigma", parameterValues);

        auto e = bd.expr();
        for( auto const& param : parameters )
        {
            if( e.expression().hasSymbol(param.first) )
            {
                if( parameters[param.first].hasMinMax() )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
                else
                    e.setParameterValues( { param.first, parameters[param.first].value() } );
            }
        }
        fV += integrate( markedfaces(M_mesh, bd.markers()),
                        e*sigma0*(M_gamma/hFace()*id(u1) - grad(u1)*N() ) );
    }

    aV.solve(_solution=u1, _rhs=fV, _name="electro");

    auto aT = form2(_test=XhT, _trial=XhT);
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto k = mat.getScalar("k", parameterValues);
        aT += integrate( markedelements(M_mesh,mat.meshMarkers()), k*inner(gradt(u2),grad(u2)) );
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : parameters )
        {
            if( e.expression().hasSymbol(param.first) )
            {
                if( parameters[param.first].hasMinMax() )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
                else
                    e.setParameterValues( { param.first, parameters[param.first].value() } );
            }
        }
        aT += integrate( markedfaces(M_mesh,bd.markers()),
                         e*inner(idt(u2),id(u2)) );
    }

    auto fT = form1(_test=XhT);
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma0 = mat.getScalar("sigma", parameterValues);
        fT += integrate( markedelements(M_mesh, mat.meshMarkers()),
                         sigma0*inner(gradv(u1))*id(u2) );
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr2();
        for( auto const& param : parameters )
        {
            if( e.expression().hasSymbol(param.first) )
            {
                if( parameters[param.first].hasMinMax() )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
                else
                    e.setParameterValues( { param.first, parameters[param.first].value() } );
            }
        }
        fT += integrate( markedfaces(M_mesh, bd.markers()),
                         e*id(u2) );
    }

    aT.solve(_solution=u2, _rhs=fT, _name="thermo");

    toc("solveLinear");

    return U;
}

double ThermoElectricNL::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    if ( need_to_solve )
        u = this->solve( mu );

    this->computeBetaQm( u, mu );


    double output=0;
    if ( output_index < Nl() )
    {
        for ( int q = 0; q < Ql(output_index); q++ )
        {
            for( int m = 0; m < mMaxF(output_index, q); ++m )
            {
                element_ptrtype eltF( new element_type( Xh ) );
                *eltF = *M_Fqm[output_index][q][m];
                output += M_betaFqm[output_index][q][m]*dot( *eltF, u );
                // output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], u );
            }
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

ThermoElectricNL::parameter_type ThermoElectricNL::parameterProperties() const
{
    auto parameters = M_modelProps->parameters();
    auto mu = Dmu->element();
    int i = 0;
    for( auto const& [key,parameter] : parameters )
        if( parameter.hasMinMax() )
            mu(i++) = parameter.value();
    return mu;
}

int ThermoElectricNL::mMaxSigma( std::string mat)
{
    auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
    return eim_sigma->mMax();
}

ThermoElectricNL::q_sigma_element_type
ThermoElectricNL::eimSigmaQ(std::string mat, int m)
{
    auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
    return eim_sigma->q(m);
}

ThermoElectricNL::vectorN_type
ThermoElectricNL::eimSigmaBeta(std::string mat, parameter_type const& mu, vectorN_type const& vtN)
{
    auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
    return eim_sigma->beta(mu, vtN);
}


}

