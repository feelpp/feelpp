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
   \file test_simplex_dynamic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-05
   \brief Tests for Simplex<Dim, Dynamic> with type-level dynamic geometric order
 */

#define BOOST_TEST_MODULE test_simplex_dynamic
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/concepts.hpp>
#include <feel/feelmesh/simplex.hpp>

#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( simplex_dynamic_suite )

//=============================================================================
// SECTION 1: Backward Compatibility Tests
// Ensure existing static-order Simplex still works exactly as before
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_static_order_backward_compat )
{
    // Test that Simplex<2, 1> works as before
    using simplex_t = Simplex<2, 1>;

    // Static assertions for compile-time properties
    static_assert( simplex_t::nDim == 2, "Expected 2D simplex" );
    static_assert( simplex_t::nOrder == 1, "Expected order 1" );
    static_assert( simplex_t::nRealDim == 2, "Expected real dim 2" );
    static_assert( simplex_t::is_simplex, "Expected is_simplex true" );
    static_assert( !simplex_t::is_hypercube, "Expected is_hypercube false" );
    static_assert( simplex_t::is_order_static, "Expected static order" );
    static_assert( !simplex_t::is_order_dynamic, "Expected not dynamic order" );

    // Static assertions for point counts (P1 triangle has 3 points)
    static_assert( simplex_t::numVertices == 3, "Triangle has 3 vertices" );
    static_assert( simplex_t::numEdges == 3, "Triangle has 3 edges" );
    static_assert( simplex_t::numPoints == 3, "P1 triangle has 3 points" );

    // Default constructor should work
    simplex_t simplex;
    BOOST_CHECK_EQUAL( simplex.topologicalDimension(), 2 );
    BOOST_CHECK_EQUAL( simplex.dimension(), 2 );
    BOOST_CHECK_EQUAL( simplex.order(), 1 );
}

BOOST_AUTO_TEST_CASE( test_simplex_static_order_p2 )
{
    // Test P2 simplex (order 2)
    using simplex_t = Simplex<2, 2>;

    static_assert( simplex_t::nOrder == 2, "Expected order 2" );
    static_assert( simplex_t::is_order_static, "Expected static order" );
    // P2 triangle: 3 vertices + 3 edge midpoints = 6 points
    static_assert( simplex_t::numPoints == 6, "P2 triangle has 6 points" );
    static_assert( simplex_t::nbPtsPerVertex == 1, "1 point per vertex" );
    static_assert( simplex_t::nbPtsPerEdge == 1, "1 point per edge for P2" );

    simplex_t simplex;
    BOOST_CHECK_EQUAL( simplex.topologicalDimension(), 2 );
    BOOST_CHECK_EQUAL( simplex.order(), 2 );
}

BOOST_AUTO_TEST_CASE( test_simplex_static_order_3d )
{
    // Test 3D simplex (tetrahedron)
    using simplex_t = Simplex<3, 1>;

    static_assert( simplex_t::nDim == 3, "Expected 3D simplex" );
    static_assert( simplex_t::numVertices == 4, "Tetra has 4 vertices" );
    static_assert( simplex_t::numEdges == 6, "Tetra has 6 edges" );
    static_assert( simplex_t::numFaces == 4, "Tetra has 4 faces" );
    static_assert( simplex_t::numPoints == 4, "P1 tetra has 4 points" );

    simplex_t simplex;
    BOOST_CHECK_EQUAL( simplex.topologicalDimension(), 3 );
}

//=============================================================================
// SECTION 2: Type-Level Dynamic RuntimeOrder Tests
// Test Simplex<Dim, Dynamic> functionality
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_type_traits )
{
    // Test that Simplex<2, Dynamic> is a distinct type with correct traits
    using simplex_static_t = Simplex<2, 1>;
    using simplex_dynamic_t = Simplex<2, Dynamic>;

    // Types should be different
    static_assert( !std::is_same_v<simplex_static_t, simplex_dynamic_t>,
                   "Static and dynamic simplex must be different types" );

    // Dynamic simplex should have is_order_dynamic = true
    static_assert( simplex_dynamic_t::is_order_dynamic,
                   "Expected is_order_dynamic true for Simplex<2, Dynamic>" );
    static_assert( !simplex_dynamic_t::is_order_static,
                   "Expected is_order_static false for Simplex<2, Dynamic>" );

    // Static simplex should have is_order_static = true
    static_assert( simplex_static_t::is_order_static,
                   "Expected is_order_static true for Simplex<2, 1>" );
    static_assert( !simplex_static_t::is_order_dynamic,
                   "Expected is_order_dynamic false for Simplex<2, 1>" );

    // Dimension should still be static
    static_assert( simplex_dynamic_t::nDim == 2, "nDim should be 2" );
    static_assert( simplex_dynamic_t::nRealDim == 2, "nRealDim should be 2" );
}

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_construction )
{
    // Test construction of dynamic-order simplex
    using simplex_dynamic_t = Simplex<2, Dynamic>;

    // Construct with RuntimeOrder(1) - P1 geometry
    simplex_dynamic_t simplex_p1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( simplex_p1.order(), 1 );
    BOOST_CHECK_EQUAL( simplex_p1.nPointsTotal(), 3 );  // P1 triangle

    // Construct with RuntimeOrder(2) - P2 geometry
    simplex_dynamic_t simplex_p2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( simplex_p2.order(), 2 );
    BOOST_CHECK_EQUAL( simplex_p2.nPointsTotal(), 6 );  // P2 triangle

    // Construct with RuntimeOrder(3) - P3 geometry
    simplex_dynamic_t simplex_p3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( simplex_p3.order(), 3 );
    BOOST_CHECK_EQUAL( simplex_p3.nPointsTotal(), 10 ); // P3 triangle
}

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_3d )
{
    // Test 3D dynamic simplex (tetrahedron)
    using simplex_dynamic_t = Simplex<3, Dynamic>;

    simplex_dynamic_t tetra_p1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( tetra_p1.order(), 1 );
    BOOST_CHECK_EQUAL( tetra_p1.nPointsTotal(), 4 );  // P1 tetra

    simplex_dynamic_t tetra_p2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( tetra_p2.order(), 2 );
    BOOST_CHECK_EQUAL( tetra_p2.nPointsTotal(), 10 ); // P2 tetra

    // Verify dimension methods still work
    BOOST_CHECK_EQUAL( tetra_p1.topologicalDimension(), 3 );
    BOOST_CHECK_EQUAL( tetra_p2.dimension(), 3 );
}

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_1d )
{
    // Test 1D dynamic simplex (segment)
    using simplex_dynamic_t = Simplex<1, Dynamic>;

    simplex_dynamic_t seg_p1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( seg_p1.order(), 1 );
    BOOST_CHECK_EQUAL( seg_p1.nPointsTotal(), 2 );  // P1 segment

    simplex_dynamic_t seg_p3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( seg_p3.order(), 3 );
    BOOST_CHECK_EQUAL( seg_p3.nPointsTotal(), 4 );  // P3 segment
}

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_0d_p0 )
{
    // Test 0D dynamic simplex (point)
    using simplex_dynamic_t = Simplex<0, Dynamic>;
    using simplex_static_t = Simplex<0, 0>;

    simplex_dynamic_t pt_p0( RuntimeOrder( 0 ) );
    simplex_static_t pt_static_p0;

    BOOST_CHECK_EQUAL( pt_p0.order(), 0 );
    BOOST_CHECK_EQUAL( pt_p0.nPointsTotal(), 1 );
    BOOST_CHECK_EQUAL( pt_p0.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( pt_p0.nPointsOnEdge(), 0 );
    BOOST_CHECK_EQUAL( pt_p0.nPointsOnFace(), 0 );
    BOOST_CHECK_EQUAL( pt_p0.nPointsOnVolume(), 0 );

    BOOST_CHECK_EQUAL( simplex_static_t::numPoints, 1 );
    BOOST_CHECK_EQUAL( simplex_static_t::nbPtsPerVertex, 1 );
    BOOST_CHECK_EQUAL( pt_static_p0.nPointsTotal(), 1 );
    BOOST_CHECK_EQUAL( pt_static_p0.nPointsOnVertex(), 1 );
}

BOOST_AUTO_TEST_CASE( test_simplex_dynamic_point_distribution )
{
    // Test per-entity point counts for dynamic order
    using simplex_dynamic_t = Simplex<2, Dynamic>;

    // P1: only vertex points
    simplex_dynamic_t p1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( p1.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( p1.nPointsOnEdge(), 0 );
    BOOST_CHECK_EQUAL( p1.nPointsOnFace(), 0 );

    // P2: vertex + edge midpoints
    simplex_dynamic_t p2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( p2.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( p2.nPointsOnEdge(), 1 );  // 1 midpoint per edge
    BOOST_CHECK_EQUAL( p2.nPointsOnFace(), 0 );

    // P3: vertex + 2 edge points + 1 face point
    simplex_dynamic_t p3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( p3.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( p3.nPointsOnEdge(), 2 );   // 2 points per edge
    BOOST_CHECK_EQUAL( p3.nPointsOnFace(), 1 );   // 1 interior face point
}

//=============================================================================
// SECTION 3: Comparison Tests - Static vs Dynamic Values
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_static_vs_dynamic_p1_2d )
{
    // Compare static Simplex<2, 1> with dynamic Simplex<2, Dynamic>(RuntimeOrder(1))
    using simplex_static_t = Simplex<2, 1>;
    using simplex_dynamic_t = Simplex<2, Dynamic>;

    simplex_static_t s_static;
    simplex_dynamic_t s_dynamic( RuntimeOrder( 1 ) );

    // RuntimeOrder should match
    BOOST_CHECK_EQUAL( s_static.order(), s_dynamic.order() );

    // Point counts should match
    BOOST_CHECK_EQUAL( s_static.nPointsTotal(), s_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnVertex(), s_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnEdge(), s_dynamic.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnFace(), s_dynamic.nPointsOnFace() );

    // Static constants should match runtime values
    BOOST_CHECK_EQUAL( simplex_static_t::numPoints, s_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( simplex_static_t::nbPtsPerVertex, s_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( simplex_static_t::nbPtsPerEdge, s_dynamic.nPointsOnEdge() );
}

BOOST_AUTO_TEST_CASE( test_simplex_static_vs_dynamic_p2_2d )
{
    // Compare static Simplex<2, 2> with dynamic Simplex<2, Dynamic>(RuntimeOrder(2))
    using simplex_static_t = Simplex<2, 2>;
    using simplex_dynamic_t = Simplex<2, Dynamic>;

    simplex_static_t s_static;
    simplex_dynamic_t s_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( s_static.order(), s_dynamic.order() );
    BOOST_CHECK_EQUAL( s_static.nPointsTotal(), s_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnVertex(), s_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnEdge(), s_dynamic.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnFace(), s_dynamic.nPointsOnFace() );
}

BOOST_AUTO_TEST_CASE( test_simplex_static_vs_dynamic_p3_2d )
{
    // Compare static Simplex<2, 3> with dynamic Simplex<2, Dynamic>(RuntimeOrder(3))
    using simplex_static_t = Simplex<2, 3>;
    using simplex_dynamic_t = Simplex<2, Dynamic>;

    simplex_static_t s_static;
    simplex_dynamic_t s_dynamic( RuntimeOrder( 3 ) );

    BOOST_CHECK_EQUAL( s_static.order(), s_dynamic.order() );
    BOOST_CHECK_EQUAL( s_static.nPointsTotal(), s_dynamic.nPointsTotal() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnVertex(), s_dynamic.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnEdge(), s_dynamic.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( s_static.nPointsOnFace(), s_dynamic.nPointsOnFace() );
}

BOOST_AUTO_TEST_CASE( test_simplex_static_vs_dynamic_3d )
{
    // Compare 3D static vs dynamic for orders 1, 2
    using simplex_static_p1_t = Simplex<3, 1>;
    using simplex_static_p2_t = Simplex<3, 2>;
    using simplex_dynamic_t = Simplex<3, Dynamic>;

    simplex_static_p1_t s_static_p1;
    simplex_static_p2_t s_static_p2;
    simplex_dynamic_t s_dynamic_p1( RuntimeOrder( 1 ) );
    simplex_dynamic_t s_dynamic_p2( RuntimeOrder( 2 ) );

    // P1 comparison
    BOOST_CHECK_EQUAL( s_static_p1.order(), s_dynamic_p1.order() );
    BOOST_CHECK_EQUAL( s_static_p1.nPointsTotal(), s_dynamic_p1.nPointsTotal() );
    BOOST_CHECK_EQUAL( s_static_p1.nPointsOnVertex(), s_dynamic_p1.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( s_static_p1.nPointsOnEdge(), s_dynamic_p1.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( s_static_p1.nPointsOnFace(), s_dynamic_p1.nPointsOnFace() );
    BOOST_CHECK_EQUAL( s_static_p1.nPointsOnVolume(), s_dynamic_p1.nPointsOnVolume() );

    // P2 comparison
    BOOST_CHECK_EQUAL( s_static_p2.order(), s_dynamic_p2.order() );
    BOOST_CHECK_EQUAL( s_static_p2.nPointsTotal(), s_dynamic_p2.nPointsTotal() );
}

//=============================================================================
// SECTION 4: RuntimeOrder Wrapper Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_order_wrapper )
{
    // Test the RuntimeOrder wrapper class
    RuntimeOrder o1( 1 );
    BOOST_CHECK_EQUAL( o1.value, 1 );

    RuntimeOrder o2( 2 );
    BOOST_CHECK_EQUAL( o2.value, 2 );

    // Test conversion to uint16_type
    uint16_type val = static_cast<uint16_type>( o1 );
    BOOST_CHECK_EQUAL( val, 1 );

    // Test construction from int
    RuntimeOrder o3( 3 );
    BOOST_CHECK_EQUAL( o3.value, 3 );
}

//=============================================================================
// SECTION 5: polyDims Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_polydims )
{
    // Test polyDims function for various dimensions and orders
    using simplex_1d_t = Simplex<1, 1>;
    using simplex_2d_t = Simplex<2, 1>;
    using simplex_3d_t = Simplex<3, 1>;

    // 1D: polyDims(n) = n + 1
    BOOST_CHECK_EQUAL( simplex_1d_t::polyDims( 0 ), 1 );
    BOOST_CHECK_EQUAL( simplex_1d_t::polyDims( 1 ), 2 );
    BOOST_CHECK_EQUAL( simplex_1d_t::polyDims( 2 ), 3 );
    BOOST_CHECK_EQUAL( simplex_1d_t::polyDims( 3 ), 4 );

    // 2D: polyDims(n) = (n+1)(n+2)/2
    BOOST_CHECK_EQUAL( simplex_2d_t::polyDims( 0 ), 1 );
    BOOST_CHECK_EQUAL( simplex_2d_t::polyDims( 1 ), 3 );
    BOOST_CHECK_EQUAL( simplex_2d_t::polyDims( 2 ), 6 );
    BOOST_CHECK_EQUAL( simplex_2d_t::polyDims( 3 ), 10 );

    // 3D: polyDims(n) = (n+1)(n+2)(n+3)/6
    BOOST_CHECK_EQUAL( simplex_3d_t::polyDims( 0 ), 1 );
    BOOST_CHECK_EQUAL( simplex_3d_t::polyDims( 1 ), 4 );
    BOOST_CHECK_EQUAL( simplex_3d_t::polyDims( 2 ), 10 );
    BOOST_CHECK_EQUAL( simplex_3d_t::polyDims( 3 ), 20 );
}

//=============================================================================
// SECTION 6: Name Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_name_static )
{
    // Test name() for static order
    using simplex_2d_p1_t = Simplex<2, 1>;
    using simplex_2d_p2_t = Simplex<2, 2>;
    using simplex_3d_p1_t = Simplex<3, 1>;

    BOOST_CHECK_EQUAL( simplex_2d_p1_t::name(), "Simplex_2_1_2" );
    BOOST_CHECK_EQUAL( simplex_2d_p2_t::name(), "Simplex_2_2_2" );
    BOOST_CHECK_EQUAL( simplex_3d_p1_t::name(), "Simplex_3_1_3" );
}

BOOST_AUTO_TEST_CASE( test_simplex_name_dynamic )
{
    // Test name() for dynamic order (instance method)
    Simplex<2, Dynamic> s1( RuntimeOrder( 1 ) );
    Simplex<2, Dynamic> s2( RuntimeOrder( 2 ) );
    Simplex<3, Dynamic> s3( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( s1.name(), "Simplex_2_1_2" );
    BOOST_CHECK_EQUAL( s2.name(), "Simplex_2_2_2" );
    BOOST_CHECK_EQUAL( s3.name(), "Simplex_3_1_3" );
}

BOOST_AUTO_TEST_CASE( test_polydims_matches_numpoints_dynamic_orders_up_to_10 )
{
    // Verify polyDims matches dynamic point count for higher runtime orders
    using simplex_2d_dyn = Simplex<2, Dynamic>;
    using simplex_3d_dyn = Simplex<3, Dynamic>;
    using simplex_1d_dyn = Simplex<1, Dynamic>;
    using simplex_2d_p1 = Simplex<2, 1>;
    using simplex_3d_p1 = Simplex<3, 1>;
    using simplex_1d_p1 = Simplex<1, 1>;

    for ( uint16_type order = 1; order <= 10; ++order )
    {
        simplex_2d_dyn s2{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( s2.nPointsTotal(), simplex_2d_p1::polyDims( order ) );

        simplex_3d_dyn s3{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( s3.nPointsTotal(), simplex_3d_p1::polyDims( order ) );

        simplex_1d_dyn s1{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( s1.nPointsTotal(), simplex_1d_p1::polyDims( order ) );
    }
}

BOOST_AUTO_TEST_SUITE_END()
