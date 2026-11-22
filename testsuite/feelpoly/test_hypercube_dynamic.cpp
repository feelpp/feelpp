/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-05

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
   \file test_hypercube_dynamic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-05
   \brief Unit tests for Hypercube<Dim, Dynamic> dynamic order support
 */

#define BOOST_TEST_MODULE test_hypercube_dynamic
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelpoly/order.hpp>

#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( hypercube_dynamic_suite )

//=============================================================================
// SECTION 1: Static Order Backward Compatibility Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_hypercube_static_2d_q1 )
{
    // Test that static Hypercube still works as expected
    using hypercube_t = Hypercube<2, 1>;

    static_assert( hypercube_t::is_order_static, "Should be static order" );
    static_assert( !hypercube_t::is_order_dynamic, "Should not be dynamic" );
    static_assert( hypercube_t::nDim == 2, "Dimension should be 2" );
    static_assert( hypercube_t::nOrder == 1, "Order should be 1" );
    static_assert( hypercube_t::numVertices == 4, "Quad has 4 vertices" );
    static_assert( hypercube_t::numEdges == 4, "Quad has 4 edges" );
    static_assert( hypercube_t::numPoints == 4, "Q1 quad has 4 points" );

    hypercube_t hypercube;
    BOOST_CHECK_EQUAL( hypercube.topologicalDimension(), 2 );
    BOOST_CHECK_EQUAL( hypercube.order(), 1 );
}

BOOST_AUTO_TEST_CASE( test_hypercube_static_2d_q2 )
{
    using hypercube_t = Hypercube<2, 2>;

    static_assert( hypercube_t::nOrder == 2, "Order should be 2" );
    // Q2: 4 vertices + 4 edge midpoints + 1 interior = 9
    static_assert( hypercube_t::numPoints == 9, "Q2 quad has 9 points" );

    hypercube_t hypercube;
    BOOST_CHECK_EQUAL( hypercube.order(), 2 );
    BOOST_CHECK_EQUAL( hypercube.nPointsTotal(), 9 );
}

BOOST_AUTO_TEST_CASE( test_hypercube_static_3d_q1 )
{
    using hypercube_t = Hypercube<3, 1>;

    static_assert( hypercube_t::nDim == 3, "Dimension should be 3" );
    static_assert( hypercube_t::numVertices == 8, "Hex has 8 vertices" );
    static_assert( hypercube_t::numEdges == 12, "Hex has 12 edges" );
    static_assert( hypercube_t::numFaces == 6, "Hex has 6 faces" );
    static_assert( hypercube_t::numPoints == 8, "Q1 hex has 8 points" );

    hypercube_t hypercube;
    BOOST_CHECK_EQUAL( hypercube.topologicalDimension(), 3 );
}

//=============================================================================
// SECTION 2: Type-Level Dynamic Order Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_type_traits )
{
    // Test that Hypercube<2, Dynamic> is a distinct type with correct traits
    using hypercube_static_t = Hypercube<2, 1>;
    using hypercube_dynamic_t = Hypercube<2, Dynamic>;

    // Type distinctness
    static_assert( !std::is_same_v<hypercube_static_t, hypercube_dynamic_t>,
                   "Static and dynamic types must be distinct" );

    // Order detection
    static_assert( hypercube_static_t::is_order_static, "Static type should be static" );
    static_assert( hypercube_dynamic_t::is_order_dynamic, "Dynamic type should be dynamic" );
    static_assert( !hypercube_dynamic_t::is_order_static, "Dynamic type should not be static" );

    // Dimension should still be compile-time
    static_assert( hypercube_dynamic_t::nDim == 2, "nDim should be 2" );
    static_assert( hypercube_dynamic_t::nRealDim == 2, "nRealDim should be 2" );
}

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_construction )
{
    // Test construction of dynamic-order hypercube
    using hypercube_dynamic_t = Hypercube<2, Dynamic>;

    // Construct with RuntimeOrder(1) - Q1 geometry
    hypercube_dynamic_t hypercube_q1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( hypercube_q1.order(), 1 );
    BOOST_CHECK_EQUAL( hypercube_q1.nPointsTotal(), 4 );  // Q1 quad: 4 vertices

    // Construct with RuntimeOrder(2) - Q2 geometry
    hypercube_dynamic_t hypercube_q2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( hypercube_q2.order(), 2 );
    BOOST_CHECK_EQUAL( hypercube_q2.nPointsTotal(), 9 );  // Q2 quad: (2+1)^2 = 9

    // Construct with RuntimeOrder(3) - Q3 geometry
    hypercube_dynamic_t hypercube_q3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( hypercube_q3.order(), 3 );
    BOOST_CHECK_EQUAL( hypercube_q3.nPointsTotal(), 16 ); // Q3 quad: (3+1)^2 = 16
}

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_3d )
{
    // Test 3D dynamic hypercube (hexahedron)
    using hypercube_dynamic_t = Hypercube<3, Dynamic>;

    hypercube_dynamic_t hex_q1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( hex_q1.order(), 1 );
    BOOST_CHECK_EQUAL( hex_q1.nPointsTotal(), 8 );  // Q1 hex: 2^3 = 8

    hypercube_dynamic_t hex_q2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( hex_q2.order(), 2 );
    BOOST_CHECK_EQUAL( hex_q2.nPointsTotal(), 27 ); // Q2 hex: 3^3 = 27

    // Verify dimension methods still work
    BOOST_CHECK_EQUAL( hex_q1.topologicalDimension(), 3 );
    BOOST_CHECK_EQUAL( hex_q2.dimension(), 3 );
}

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_1d )
{
    // Test 1D dynamic hypercube (segment)
    using hypercube_dynamic_t = Hypercube<1, Dynamic>;

    hypercube_dynamic_t seg_q1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( seg_q1.order(), 1 );
    BOOST_CHECK_EQUAL( seg_q1.nPointsTotal(), 2 );  // Q1 segment: 2 points

    hypercube_dynamic_t seg_q3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( seg_q3.order(), 3 );
    BOOST_CHECK_EQUAL( seg_q3.nPointsTotal(), 4 );  // Q3 segment: 4 points
}

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_0d_q0 )
{
    // Test 0D dynamic hypercube (point)
    using hypercube_dynamic_t = Hypercube<0, Dynamic>;
    using hypercube_static_t = Hypercube<0, 0>;

    hypercube_dynamic_t pt_q0( RuntimeOrder( 0 ) );
    hypercube_static_t pt_static_q0;

    BOOST_CHECK_EQUAL( pt_q0.order(), 0 );
    BOOST_CHECK_EQUAL( pt_q0.nPointsTotal(), 1 );
    BOOST_CHECK_EQUAL( pt_q0.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( pt_q0.nPointsOnEdge(), 0 );
    BOOST_CHECK_EQUAL( pt_q0.nPointsOnFace(), 0 );
    BOOST_CHECK_EQUAL( pt_q0.nPointsOnVolume(), 0 );

    BOOST_CHECK_EQUAL( hypercube_static_t::numPoints, 1 );
    BOOST_CHECK_EQUAL( hypercube_static_t::nbPtsPerVertex, 1 );
    BOOST_CHECK_EQUAL( pt_static_q0.nPointsTotal(), 1 );
    BOOST_CHECK_EQUAL( pt_static_q0.nPointsOnVertex(), 1 );
}

BOOST_AUTO_TEST_CASE( test_hypercube_dynamic_point_distribution )
{
    // Test per-entity point counts for dynamic order
    using hypercube_dynamic_t = Hypercube<2, Dynamic>;

    // Q1: only vertex points
    hypercube_dynamic_t q1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( q1.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( q1.nPointsOnEdge(), 0 );
    BOOST_CHECK_EQUAL( q1.nPointsOnFace(), 0 );

    // Q2: vertex + edge midpoints + 1 face point
    hypercube_dynamic_t q2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( q2.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( q2.nPointsOnEdge(), 1 );   // 1 midpoint per edge
    BOOST_CHECK_EQUAL( q2.nPointsOnFace(), 1 );   // 1 interior face point

    // Q3: vertex + 2 edge points + 4 face points
    hypercube_dynamic_t q3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( q3.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( q3.nPointsOnEdge(), 2 );   // 2 points per edge
    BOOST_CHECK_EQUAL( q3.nPointsOnFace(), 4 );   // (3-1)^2 = 4 interior face points
}

//=============================================================================
// SECTION 3: Comparison Tests - Static vs Dynamic Values
//=============================================================================

BOOST_AUTO_TEST_CASE( test_hypercube_static_vs_dynamic_q1_2d )
{
    // Compare static Hypercube<2, 1> with dynamic Hypercube<2, Dynamic>(RuntimeOrder(1))
    using hypercube_static_t = Hypercube<2, 1>;
    using hypercube_dynamic_t = Hypercube<2, Dynamic>;

    hypercube_static_t h_static;
    hypercube_dynamic_t h_dynamic( RuntimeOrder( 1 ) );

    // Order should match
    BOOST_CHECK_EQUAL( h_static.order(), h_dynamic.order() );

    // Point counts should match
    BOOST_CHECK_EQUAL( h_static.nPointsTotal(), h_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnVertex(), h_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnEdge(), h_dynamic.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnFace(), h_dynamic.nPointsOnFace() );

    // Static constants should match runtime values
    BOOST_CHECK_EQUAL( hypercube_static_t::numPoints, h_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( hypercube_static_t::nbPtsPerVertex, h_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( hypercube_static_t::nbPtsPerEdge, h_dynamic.nPointsOnEdge() );
}

BOOST_AUTO_TEST_CASE( test_hypercube_static_vs_dynamic_q2_2d )
{
    // Compare static Hypercube<2, 2> with dynamic Hypercube<2, Dynamic>(RuntimeOrder(2))
    using hypercube_static_t = Hypercube<2, 2>;
    using hypercube_dynamic_t = Hypercube<2, Dynamic>;

    hypercube_static_t h_static;
    hypercube_dynamic_t h_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( h_static.order(), h_dynamic.order() );
    BOOST_CHECK_EQUAL( h_static.nPointsTotal(), h_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnVertex(), h_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnEdge(), h_dynamic.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnFace(), h_dynamic.nPointsOnFace() );
}

BOOST_AUTO_TEST_CASE( test_hypercube_static_vs_dynamic_q1_3d )
{
    // Compare 3D static vs dynamic for Q1
    using hypercube_static_t = Hypercube<3, 1>;
    using hypercube_dynamic_t = Hypercube<3, Dynamic>;

    hypercube_static_t h_static;
    hypercube_dynamic_t h_dynamic( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( h_static.order(), h_dynamic.order() );
    BOOST_CHECK_EQUAL( h_static.nPointsTotal(), h_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnVertex(), h_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnEdge(), h_dynamic.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnFace(), h_dynamic.nPointsOnFace() );
    BOOST_CHECK_EQUAL( h_static.nPointsOnVolume(), h_dynamic.nPointsOnVolume() );
}

BOOST_AUTO_TEST_CASE( test_hypercube_static_vs_dynamic_q2_3d )
{
    // Compare 3D static vs dynamic for Q2
    using hypercube_static_t = Hypercube<3, 2>;
    using hypercube_dynamic_t = Hypercube<3, Dynamic>;

    hypercube_static_t h_static;
    hypercube_dynamic_t h_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( h_static.order(), h_dynamic.order() );
    BOOST_CHECK_EQUAL( h_static.nPointsTotal(), h_dynamic.nPointsTotal() );
}

//=============================================================================
// SECTION 4: polyDims Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_hypercube_polydims )
{
    // Test polyDims function for various dimensions and orders
    using hypercube_1d_t = Hypercube<1, 1>;
    using hypercube_2d_t = Hypercube<2, 1>;
    using hypercube_3d_t = Hypercube<3, 1>;

    // 1D: polyDims(n) = n + 1
    BOOST_CHECK_EQUAL( hypercube_1d_t::polyDims( 0 ), 1 );
    BOOST_CHECK_EQUAL( hypercube_1d_t::polyDims( 1 ), 2 );
    BOOST_CHECK_EQUAL( hypercube_1d_t::polyDims( 2 ), 3 );
    BOOST_CHECK_EQUAL( hypercube_1d_t::polyDims( 3 ), 4 );

    // 2D: polyDims(n) = (n+1)^2
    BOOST_CHECK_EQUAL( hypercube_2d_t::polyDims( 0 ), 1 );
    BOOST_CHECK_EQUAL( hypercube_2d_t::polyDims( 1 ), 4 );
    BOOST_CHECK_EQUAL( hypercube_2d_t::polyDims( 2 ), 9 );
    BOOST_CHECK_EQUAL( hypercube_2d_t::polyDims( 3 ), 16 );

    // 3D: polyDims(n) = (n+1)^3
    BOOST_CHECK_EQUAL( hypercube_3d_t::polyDims( 0 ), 1 );
    BOOST_CHECK_EQUAL( hypercube_3d_t::polyDims( 1 ), 8 );
    BOOST_CHECK_EQUAL( hypercube_3d_t::polyDims( 2 ), 27 );
    BOOST_CHECK_EQUAL( hypercube_3d_t::polyDims( 3 ), 64 );
}

//=============================================================================
// SECTION 5: Name Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_hypercube_name_static )
{
    // Test name() for static order
    using hypercube_2d_q1_t = Hypercube<2, 1>;
    using hypercube_2d_q2_t = Hypercube<2, 2>;
    using hypercube_3d_q1_t = Hypercube<3, 1>;

    BOOST_CHECK_EQUAL( hypercube_2d_q1_t::name(), "Hypercube_2_1_2" );
    BOOST_CHECK_EQUAL( hypercube_2d_q2_t::name(), "Hypercube_2_2_2" );
    BOOST_CHECK_EQUAL( hypercube_3d_q1_t::name(), "Hypercube_3_1_3" );
}

BOOST_AUTO_TEST_CASE( test_hypercube_name_dynamic )
{
    // Test name() for dynamic order (instance method)
    Hypercube<2, Dynamic> h1( RuntimeOrder( 1 ) );
    Hypercube<2, Dynamic> h2( RuntimeOrder( 2 ) );
    Hypercube<3, Dynamic> h3( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( h1.name(), "Hypercube_2_1_2" );
    BOOST_CHECK_EQUAL( h2.name(), "Hypercube_2_2_2" );
    BOOST_CHECK_EQUAL( h3.name(), "Hypercube_3_1_3" );
}

//=============================================================================
// SECTION 6: Cross-validation with polyDims
//=============================================================================

BOOST_AUTO_TEST_CASE( test_polydims_matches_numpoints )
{
    // Verify polyDims matches numPoints for various orders
    using hypercube_2d_dyn = Hypercube<2, Dynamic>;
    using hypercube_3d_dyn = Hypercube<3, Dynamic>;
    using hypercube_1d_dyn = Hypercube<1, Dynamic>;
    using hypercube_2d_q1 = Hypercube<2, 1>;
    using hypercube_3d_q1 = Hypercube<3, 1>;
    using hypercube_1d_q1 = Hypercube<1, 1>;

    for ( uint16_type order = 1; order <= 10; ++order )
    {
        // 2D
        hypercube_2d_dyn h2{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( h2.nPointsTotal(), hypercube_2d_q1::polyDims( order ) );

        // 3D
        hypercube_3d_dyn h3{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( h3.nPointsTotal(), hypercube_3d_q1::polyDims( order ) );

        // 1D
        hypercube_1d_dyn h1{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( h1.nPointsTotal(), hypercube_1d_q1::polyDims( order ) );
    }
}

BOOST_AUTO_TEST_SUITE_END()
