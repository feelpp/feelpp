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
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::sparse_matrix_ptrtype
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::assembleForMDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto mesh = this->Xh->mesh();
    auto u1 = u.template element<0>();
    auto u2 = u.template element<1>();
    auto uOld1 = u.template element<0>();
    auto uOld2 = u.template element<1>();

    auto a = form2(_test=this->Xh, _trial=this->Xh);
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
        auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
        auto k = sigma*mu.parameterNamed("L")*idv(uOld2);
        a += integrate( markedelements(mesh, mat.meshMarkers()),
                        k*inner(gradt(u2),grad(u2)) );
    }
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
        auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
        a += integrate( markedelements(mesh, mat.meshMarkers()),
                            sigma*inner(gradt(u1),grad(u1)) );
    }

    // Dirichlet condition
    for( auto const&[bdName, bd] : potentialDirichlet )
    {
        auto mat = M_elecMaterials[bd.material()];
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
        auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
        a += integrate( markedfaces(mesh,bd.markers()),
                        sigma*(M_gamma/hFace()*inner(idt(u1),id(u1))
                               -inner(gradt(u1)*N(),id(u1))
                               -inner(grad(u1)*N(),idt(u1)) ) );
    }

    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr1();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        a += integrate( markedfaces(mesh,bd.markers()),
                        e.evaluate()(0,0)*inner(idt(u2),id(u2)) );
    }

    auto am = a.matrixPtr();
    am->close();

    return am;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::vector_ptrtype
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto mesh = this->Xh->mesh();
    auto u1 = u.template element<0>();
    auto u2 = u.template element<1>();
    auto uOld1 = u.template element<0>();
    auto uOld2 = u.template element<1>();

    auto f = form1(_test=this->Xh);
    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto sigma0 = mu.parameterNamed(mat.getString("misc.sigmaKey"));
        auto alpha = mu.parameterNamed(mat.getString("misc.alphaKey"));
        auto sigma = sigma0/(cst(1.)+alpha*(idv(uOld2)-cst(293.)));
        f += integrate( markedelements(mesh,mat.meshMarkers()),
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

        f += integrate( markedfaces(mesh, bd.markers()),
                        e.evaluate()(0,0)*sigma*(M_gamma/hFace()*id(u1) - grad(u1)*N() ) );
    }
    for( auto const&[bdName, bd] : temperatureRobin )
    {
        auto e = bd.expr2();
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        f += integrate( markedfaces(mesh, bd.markers()),
                        e.evaluate()(0,0)*id(u2) );
    }

    auto fv = f.vectorPtr();
    fv->close();

    return fv;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::element_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::solve(parameter_type const& mu)
{
    tic();

    auto U = this->Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();
    auto UOld = this->Xh->element();
    auto uOld1 = UOld.template element<0>();
    auto uOld2 = UOld.template element<1>();
    uOld1 = M_InitialGuess[0][0]->template element<0>();
    uOld2 = M_InitialGuess[1][0]->template element<1>();

    double increment = 0;
    int it = 0;

    do
    {
        auto am = this->assembleForMDEIMnl(mu, UOld, 0);
        auto a = form2(_test=this->Xh, _trial=this->Xh, _matrix=am);
        auto fv = this->assembleForDEIMnl(mu, UOld, 0);
        auto f = form1(_test=this->Xh, _vector=fv);
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

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::element_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::solveLinear(parameter_type const& mu)
{
    tic();
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");

    auto XhV = this->Xh->template functionSpace<0>();
    auto XhT = this->Xh->template functionSpace<1>();
    auto U = this->Xh->element();
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
                         sigma0*(M_gamma/hFace()*inner(idt(u1),id(u1))
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
                        e.evaluate()(0,0)*sigma0*(M_gamma/hFace()*id(u1) - grad(u1)*N() ) );
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
        auto e = bd.expr2();
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

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
double
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
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
                element_ptrtype eltF( new element_type( this->Xh ) );
                *eltF = *this->M_Fqm[output_index][q][m];
                output += this->M_betaFqm[output_index][q][m]*dot( *eltF, u );
                // output += this->M_betaFqm[output_index][q][m]*dot( this->M_Fqm[output_index][q][m], u );
            }
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::parameter_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::parameterProperties() const
{
    auto parameters = M_modelProps->parameters();
    auto mu = this->Dmu->element();
    int i = 0;
    for( auto const& [key,parameter] : parameters )
        if( parameter.hasMinMax() )
            mu(i++) = parameter.value();
    return mu;
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
int
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::mMaxSigma( std::string mat)
{
    auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
    return eim_sigma->mMax();
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::q_sigma_element_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::eimSigmaQ(std::string mat, int m)
{
    auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
    return eim_sigma->q(m);
}

THERMOELECTRICNL_CLASS_TEMPLATE_DECLARATIONS
typename THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::vectorN_type
THERMOELECTRICNL_CLASS_TEMPLATE_TYPE::eimSigmaBeta(std::string mat, parameter_type const& mu, vectorN_type const& vtN)
{
    auto eim_sigma = this->scalarContinuousEim()[M_elecEimIndex[mat]];
    return eim_sigma->beta(mu, vtN);
}


}

