/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-02-17

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
 * @file test_dynamic_order_math_parity.cpp
 * @brief Dynamic-order mathematical parity checks (static vs dynamic)
 */
#define BOOST_TEST_MODULE test_dynamic_order_math_parity
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdhm.hpp>
#include <feel/feeldiscr/pdhv.hpp>
#include <feel/feelvf/vf.hpp>

#include <string>
#include <type_traits>
#include <utility>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using mesh_type = Mesh<Simplex<2,1>>;
using mesh_ptrtype = std::shared_ptr<mesh_type>;

mesh_ptrtype
makeOneElementSimplexMesh( std::string const& name, bool referenceDomain )
{
    auto const domainType = referenceDomain ? Gmsh::GMSH_REFERENCE_DOMAIN : Gmsh::GMSH_REAL_DOMAIN;
    auto desc = std::make_shared<GmshSimplexDomain>( 2, 1, domainType );
    desc->setPrefix( name );
    desc->usePhysicalNames( true );
    if ( referenceDomain )
    {
        desc->setCharacteristicLength( 2.0 );
    }
    else
    {
        desc->setX( { 0.2, 1.4 } );
        desc->setY( { -0.1, 1.1 } );
        desc->setCharacteristicLength( 1.2 );
    }
    return createGMSHMesh( _mesh = new mesh_type, _desc = desc, _h = desc->h() );
}

template<int Order>
auto manufacturedScalarExpr()
{
    if constexpr ( Order == 0 )
        return cst( 1.0 );
    else if constexpr ( Order == 1 )
        return cst( 1.0 ) + 2.0 * Px() - 0.5 * Py();
    else if constexpr ( Order == 2 )
        return cst( 1.0 ) + Px() + Py() + Px() * Px() + Px() * Py() + Py() * Py();
    else if constexpr ( Order == 3 )
        return cst( 1.0 ) + Px() * Px() * Px() + Px() * Px() * Py() + Px() * Py() * Py() + Py() * Py() * Py();
    else
        return cst( 1.0 ) + pow( Px(), 4 ) + pow( Py(), 4 ) + Px() * Px() * Py() * Py();
}

template<int Order>
auto manufacturedVectorExpr()
{
    if constexpr ( Order == 0 )
        return vec( cst( 0.75 ), cst( -1.25 ) );
    else if constexpr ( Order == 1 )
        return vec( cst( 0.75 ) + 2.0 * Px() - Py(),
                    cst( -1.25 ) + Px() + 3.0 * Py() );
    else if constexpr ( Order == 2 )
        return vec( cst( 0.75 ) + Px() * Px() + Px() * Py(),
                    cst( -1.25 ) + Py() * Py() + 2.0 * Px() * Py() );
    else if constexpr ( Order == 3 )
        return vec( cst( 0.75 ) + pow( Px(), 3 ) + Px() * Py() * Py(),
                    cst( -1.25 ) + pow( Py(), 3 ) + Px() * Px() * Py() );
    else
        return vec( cst( 0.75 ) + pow( Px(), 4 ) + Px() * Py() * Py() * Py(),
                    cst( -1.25 ) + pow( Py(), 4 ) + Px() * Px() * Px() * Py() );
}

template<int Order>
auto manufacturedTensorExpr()
{
    if constexpr ( Order == 0 )
        return mat<2,2>( cst( 1.0 ), cst( -2.0 ), cst( 0.5 ), cst( 3.0 ) );
    else if constexpr ( Order == 1 )
        return mat<2,2>( cst( 1.0 ) + Px(),
                         cst( -2.0 ) + Py(),
                         cst( 0.5 ) + Px() - Py(),
                         cst( 3.0 ) + Px() + Py() );
    else if constexpr ( Order == 2 )
        return mat<2,2>( cst( 1.0 ) + Px() * Px() + Px() * Py(),
                         cst( -2.0 ) + Py() * Py(),
                         cst( 0.5 ) + Px() * Py(),
                         cst( 3.0 ) + Px() * Px() + Py() * Py() );
    else if constexpr ( Order == 3 )
        return mat<2,2>( cst( 1.0 ) + pow( Px(), 3 ),
                         cst( -2.0 ) + Px() * Px() * Py(),
                         cst( 0.5 ) + Px() * Py() * Py(),
                         cst( 3.0 ) + pow( Py(), 3 ) );
    else
        return mat<2,2>( cst( 1.0 ) + pow( Px(), 4 ),
                         cst( -2.0 ) + Px() * Px() * Py() * Py(),
                         cst( 0.5 ) + Px() * pow( Py(), 3 ),
                         cst( 3.0 ) + pow( Py(), 4 ) );
}

template<int Order>
void
checkScalarParity( mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pch<Order>( mesh );
    auto VhDynamic = Pch<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto expr = manufacturedScalarExpr<Order>();

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_REQUIRE_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto uStatic = project( _space = VhStatic, _expr = expr );
    auto uDynamic = project( _space = VhDynamic, _expr = expr );

    auto const errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - expr, _quad = _Q<20>() );
    auto const errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - expr, _quad = _Q<20>() );
    auto const errParity = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<20>() );

    BOOST_TEST_CONTEXT( "Pch static/dynamic parity P" << Order )
    {
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errParity, tol );
    }
}

template<int Order>
void
checkVectorParity( mesh_ptrtype const& mesh, double tol )
{
    auto VhStatic = Pdhv<Order>( mesh );
    auto VhDynamic = Pdhv<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto expr = manufacturedVectorExpr<Order>();

    BOOST_REQUIRE_EQUAL( VhDynamic->runtimeOrder(), Order );
    BOOST_REQUIRE_EQUAL( VhStatic->nDof(), VhDynamic->nDof() );

    auto uStatic = VhStatic->element();
    auto uDynamic = VhDynamic->element();
    uStatic.on( _range = elements( mesh ), _expr = expr );
    uDynamic.on( _range = elements( mesh ), _expr = expr );

    auto const errStatic = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - expr, _quad = _Q<20>() );
    auto const errDynamic = normL2( _range = elements( mesh ), _expr = idv( uDynamic ) - expr, _quad = _Q<20>() );
    auto const errParity = normL2( _range = elements( mesh ), _expr = idv( uStatic ) - idv( uDynamic ), _quad = _Q<20>() );

    BOOST_TEST_CONTEXT( "Pdhv static/dynamic parity P" << Order )
    {
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errParity, tol );
    }
}

template<int Order>
void
checkTensorParity( mesh_ptrtype const& mesh, double tol )
{
    auto MhStatic = Pdhm<Order>( mesh );
    auto MhDynamic = Pdhm<Dynamic>( mesh, RuntimeOrder{ Order } );
    auto expr = manufacturedTensorExpr<Order>();

    BOOST_REQUIRE_EQUAL( MhDynamic->runtimeOrder(), Order );
    BOOST_REQUIRE_EQUAL( MhStatic->nDof(), MhDynamic->nDof() );

    auto mStatic = MhStatic->element();
    auto mDynamic = MhDynamic->element();
    mStatic.on( _range = elements( mesh ), _expr = expr );
    mDynamic.on( _range = elements( mesh ), _expr = expr );

    auto const errStatic = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - expr, _quad = _Q<20>() );
    auto const errDynamic = normL2( _range = elements( mesh ), _expr = idv( mDynamic ) - expr, _quad = _Q<20>() );
    auto const errParity = normL2( _range = elements( mesh ), _expr = idv( mStatic ) - idv( mDynamic ), _quad = _Q<20>() );

    BOOST_TEST_CONTEXT( "Pdhm static/dynamic parity P" << Order )
    {
        BOOST_CHECK_SMALL( errStatic, tol );
        BOOST_CHECK_SMALL( errDynamic, tol );
        BOOST_CHECK_SMALL( errParity, tol );
    }
}

template<int... Orders, typename Fn>
void
forEachOrder( std::integer_sequence<int, Orders...>, Fn&& fn )
{
    ( fn( std::integral_constant<int, Orders>{} ), ... );
}

template<int... Orders>
void
runDynamicOrderChecks( mesh_ptrtype const& mesh )
{
    constexpr double tolScalar = 1e-10;
    constexpr double tolVector = 1e-10;
    constexpr double tolTensor = 1e-10;

    forEachOrder( std::integer_sequence<int, Orders...>{},
                  [&]( auto orderC )
                  {
                      constexpr int order = decltype( orderC )::value;
                      checkScalarParity<order>( mesh, tolScalar );
                      checkVectorParity<order>( mesh, tolVector );
                      checkTensorParity<order>( mesh, tolTensor );
                  } );
}
} // namespace

BOOST_AUTO_TEST_SUITE( dynamic_order_math_parity_suite )

BOOST_AUTO_TEST_CASE( dynamic_order_reference_simplex_math_parity )
{
    auto mesh = makeOneElementSimplexMesh( "dynamic_order_reference", true );
    runDynamicOrderChecks<0,1,2,3,4>( mesh );
}

BOOST_AUTO_TEST_CASE( dynamic_order_affine_simplex_math_parity )
{
    auto mesh = makeOneElementSimplexMesh( "dynamic_order_affine", false );
    runDynamicOrderChecks<0,1,2,3,4>( mesh );
}

BOOST_AUTO_TEST_SUITE_END()
