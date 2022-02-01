/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 08 Jan 2020

 Copyright (C) 2020 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_CG_LAPLACIAN_HPP
#define FEELPP_CG_LAPLACIAN_HPP 1

#include <feel/feelmodels/modelproperties.hpp>

namespace Feel {

/**
 * Basic laplacian implementation with Dirichlet condition
 * \tparam FunctionSpace type
 * \tparam DataT data type for right hand side, source and diffusion terms
 * 
 * The following code solve \f$-\Delta u = 1, u = 0 \mbox{ on } \Omega \f$ using
 * \f$P1\f$ piecewise polynomial continuous functions
 * \code
 * cgLaplacianDirichlet( Pch<1>(mesh), std::tuple{ expr("1"), expr("1"),expr("0") } );
 * \endcode
 */
template<typename SpaceType, typename DataT, typename = std::enable_if_t<is_functionspace_v<SpaceType>>>
auto
cgLaplacianDirichlet( std::shared_ptr<SpaceType> const& Vh, DataT && data  )
{
    auto [k,f,g] = data;
    // tag::v[]
    auto u = Vh->element();
    auto v = Vh->element( g, "g" );
    // end::v[]

    // tag::forms[]
    auto l = form1( _test = Vh );
    l = integrate( _range = elements( support( Vh ) ),
                   _expr = f * id( v ) );

    auto a = form2( _trial = Vh, _test = Vh );
    a = integrate( _range = elements( support( Vh ) ),
                   _expr = inner( k*gradt( u ), grad( v ) ) );
    a += on( _range = boundaryfaces( support( Vh ) ), _rhs = l, _element = u, _expr = g );
    // end::forms[]

    // tag::solve[]
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    a.solve( _rhs = l, _solution = u );
    // end::solve[]
    return u;
}

template<typename SpaceType, typename DataT, typename = std::enable_if_t<is_functionspace_v<SpaceType>>>
auto
cgLaplacian( std::shared_ptr<SpaceType> const& Vh, DataT && data  )
{
    auto [k,f,g,un,r_1,r_2] = data;
    tic();
    // tag::v[]
    auto u = Vh->element();
    auto v = Vh->element( g, "g" );
    // end::v[]
    toc( "Vh elements" );
    // end::mesh_space[]

    // tag::forms[]
    tic();
    auto l = form1( _test = Vh );
    l = integrate( _range = elements( support( Vh ) ),
                   _expr = f * id( v ) );
    l += integrate( _range = markedfaces( support( Vh ), "Robin" ), _expr = -r_2 * id( v ) );
    l += integrate( _range = markedfaces( support( Vh ), "Neumann" ), _expr = -un * id( v ) );
    toc( "l" );

    tic();
    auto a = form2( _trial = Vh, _test = Vh );
    tic();
    a = integrate( _range = elements( support( Vh ) ),
                   _expr = inner( k * gradt( u ), grad( v ) ) );
    toc( "a.gradgrad" );
    tic();
    a += integrate( _range = markedfaces( support( Vh ), "Robin" ), _expr = r_1 * idt( u ) * id( v ), _quad = 4 );
    toc( "a.robin" );
    tic();
    a += on( _range = markedfaces( support( Vh ), "Dirichlet" ), _rhs = l, _element = u, _expr = g );
    //! if no markers Robin Neumann or Dirichlet are present in the mesh then
    //! impose Dirichlet boundary conditions over the entire boundary
    if ( !support( Vh )->hasAnyMarker( {"Robin", "Neumann", "Dirichlet"} ) )
        a += on( _range = boundaryfaces( support( Vh ) ), _rhs = l, _element = u, _expr = g );
    toc( "a.dirichlet" );
    toc( "a" );
    // end::forms[]

    // tag::solve[]
    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve( _rhs = l, _solution = u );
    toc( "a.solve" );
    // end::solve[]
    return boption( "no-solve" )?std::nullopt:std::optional{u};
}

template<typename SpaceType, typename = std::enable_if_t<is_functionspace_v<SpaceType>>>
auto
cgLaplacianModel( std::shared_ptr<SpaceType> const& Vh, ModelProperties const& props )
{
    tic();
    // tag::v[]
    auto u = Vh->element();
    auto v = Vh->element();
    // end::v[]
    toc( "Vh elements" );
    // end::mesh_space[]

    // tag::forms[]
    tic();
    auto l = form1( _test = Vh );

    for ( auto const& [ bcid, bc ] : props.boundaryConditions2().flatten() )
    {
        
        if ( bcid.type() == "Load" )
        {
            l += integrate( _range = markedelements( support( Vh ), bc.markers() ), 
                                _expr = expr(bc.expression()) * id( v ) );
        }
        if ( bcid.type() == "Neumann" )
        {
            l += integrate( _range = markedfaces( support( Vh ), bc.markers() ), 
                            _expr = -expr(bc.expression()) * id( v ) );
        }
        if ( bcid.type() == "Robin" )
        {
            l += integrate( _range = markedfaces( support( Vh ), bc.markers() ), 
                            _expr = -expr(bc.expression2()) * id( v ) );
        }
    }
    toc( "l" );

    tic();
    auto a = form2( _trial = Vh, _test = Vh );
    tic();
    for ( auto const& [ name, mat ] : props.materials() )
    {
        if ( mat.hasProperty( "k" ) )
        {
            Feel::cout << " - Material " << name << std::endl;
            auto k = mat.property( "k" ).template expr<1,1>();
            a = integrate( _range = markedelements( support( Vh ), mat.meshMarkers() ),
                           _expr = inner( k * gradt( u ), grad( v ) ) );
        }
        else
            throw std::logic_error( std::string("missing diffusion coefficient k in material ") + name );
    }
    toc( "a.gradgrad" );
    for ( auto const& [ bcid, bc ] : props.boundaryConditions2().flatten() )
    {
        if ( bcid.type() == "Robin" )
        {
            Feel::cout << " - Boundary Condition " << bcid << std::endl;
            a += integrate( _range = markedfaces( support( Vh ), bc.markers() ), 
                            _expr = expr(bc.expression1()) * idt(v) * id( v ) );
        }
    }
    // do Dirichlet conditions last
    for ( auto const& [ bcid, bc ] : props.boundaryConditions2().flatten() )
    {
        if ( bcid.type() == "Dirichlet" )
            a += on( _range = markedfaces( support( Vh ), bc.markers() ), _rhs = l, _element = u, _expr = expr(bc.expression())  );
    }
    toc( "a" );
    // end::forms[]

    // tag::solve[]
    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    a.solve( _rhs = l, _solution = u );
    toc( "a.solve" );

    auto to_export = Vh->elementsMap();
    to_export["u"] = u;

    // end::solve[]
    return to_export;
}


}
#endif
