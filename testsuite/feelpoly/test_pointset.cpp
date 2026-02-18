/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2026-01-06

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
   \file test_pointset.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-06
   \brief Tests for PointSetEquiSpaced with static and dynamic order
 */

#define BOOST_TEST_MODULE test_pointset
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelpoly/order.hpp>

#include <cmath>
#include <string>
#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( pointset_equispaced_suite )

//=============================================================================
// Helper function to compare point matrices
//=============================================================================

template<typename MatrixType>
bool comparePointMatrices( const MatrixType& m1, const MatrixType& m2, double tol = 1e-12 )
{
    if ( m1.size1() != m2.size1() || m1.size2() != m2.size2() )
        return false;

    for ( size_t i = 0; i < m1.size1(); ++i )
    {
        for ( size_t j = 0; j < m1.size2(); ++j )
        {
            if ( std::abs( m1( i, j ) - m2( i, j ) ) > tol )
                return false;
        }
    }
    return true;
}

template<int Order, typename ConvexType>
void checkStaticDynamicPointsetEquality( double tol = 1e-12 )
{
    using pointset_static_t = PointSetEquiSpaced<ConvexType, Order, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<ConvexType, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic{ RuntimeOrder{ static_cast<uint16_type>( Order ) } };

    BOOST_TEST_CONTEXT( "static/dynamic pointset equality order=" << Order )
    {
        BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
        BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
        BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points(), tol ) );
    }
}

template<int MaxOrder, typename ConvexType, int CurrentOrder = 0>
void checkStaticDynamicPointsetRange( double tol = 1e-12 )
{
    checkStaticDynamicPointsetEquality<CurrentOrder, ConvexType>( tol );
    if constexpr ( CurrentOrder < MaxOrder )
        checkStaticDynamicPointsetRange<MaxOrder, ConvexType, CurrentOrder+1>( tol );
}

//=============================================================================
// SECTION 1: Static Order Tests (Backward Compatibility)
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_static_p1_1d )
{
    // Test P1 segment (1D)
    using convex_t = Simplex<1, 1, 1>;
    using pointset_t = PointSetEquiSpaced<convex_t, 1, double>;

    static_assert( pointset_t::is_order_static, "Expected static order" );
    static_assert( !pointset_t::is_order_dynamic, "Expected not dynamic" );
    static_assert( pointset_t::numPoints == 2, "P1 segment has 2 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 2 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 2 );
}

BOOST_AUTO_TEST_CASE( test_pointset_construction_traits )
{
    using convex_t = Simplex<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 2, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    static_assert( pointset_static_t::is_order_static, "Expected static pointset order" );
    static_assert( pointset_dynamic_t::is_order_dynamic, "Expected dynamic pointset order" );
    static_assert( std::is_constructible_v<pointset_dynamic_t, RuntimeOrder>,
                   "Dynamic pointset must be constructible from RuntimeOrder" );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p2_1d )
{
    // Test P2 segment (1D)
    using convex_t = Simplex<1, 1, 1>;
    using pointset_t = PointSetEquiSpaced<convex_t, 2, double>;

    static_assert( pointset_t::numPoints == 3, "P2 segment has 3 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p1_2d )
{
    // Test P1 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 1, double>;

    static_assert( pointset_t::numPoints == 3, "P1 triangle has 3 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size1(), 2 );  // 2D coordinates
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );  // 3 points
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p2_2d )
{
    // Test P2 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 2, double>;

    static_assert( pointset_t::numPoints == 6, "P2 triangle has 6 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 6 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 6 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p3_2d )
{
    // Test P3 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 3, double>;

    static_assert( pointset_t::numPoints == 10, "P3 triangle has 10 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 3 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 10 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p1_3d )
{
    // Test P1 tetrahedron (3D)
    using convex_t = Simplex<3, 1, 3>;
    using pointset_t = PointSetEquiSpaced<convex_t, 1, double>;

    static_assert( pointset_t::numPoints == 4, "P1 tetra has 4 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 4 );
    BOOST_CHECK_EQUAL( pts.points().size1(), 3 );  // 3D coordinates
    BOOST_CHECK_EQUAL( pts.points().size2(), 4 );  // 4 points
}

BOOST_AUTO_TEST_CASE( test_pointset_static_q1_2d )
{
    // Test Q1 quadrilateral (2D hypercube)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 1, double>;

    static_assert( pointset_t::numPoints == 4, "Q1 quad has 4 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 4 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_q2_2d )
{
    // Test Q2 quadrilateral (2D hypercube)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 2, double>;

    static_assert( pointset_t::numPoints == 9, "Q2 quad has 9 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 9 );
}

//=============================================================================
// SECTION 2: Dynamic Order Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p1_1d )
{
    // Test dynamic P1 segment (1D)
    using convex_t = Simplex<1, 1, 1>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    static_assert( pointset_t::is_order_dynamic, "Expected dynamic order" );
    static_assert( !pointset_t::is_order_static, "Expected not static" );

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 2 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 2 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p2_1d )
{
    // Test dynamic P2 segment (1D)
    using convex_t = Simplex<1, 1, 1>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p1_2d )
{
    // Test dynamic P1 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size1(), 2 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p2_2d )
{
    // Test dynamic P2 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 6 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p3_2d )
{
    // Test dynamic P3 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 3 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 10 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p1_3d )
{
    // Test dynamic P1 tetrahedron (3D)
    using convex_t = Simplex<3, 1, 3>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 4 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_q1_2d )
{
    // Test dynamic Q1 quadrilateral (2D)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 4 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_q2_2d )
{
    // Test dynamic Q2 quadrilateral (2D)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 9 );
}

//=============================================================================
// SECTION 3: Static vs Dynamic Comparison - Point Matrix Equality
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_compare_p1_1d )
{
    // Compare static and dynamic P1 segment - points must be identical
    using convex_t = Simplex<1, 1, 1>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 1, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 1 ) );

    // Check dimensions match
    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );

    // Check point matrices are identical
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p2_1d )
{
    // Compare static and dynamic P2 segment
    using convex_t = Simplex<1, 1, 1>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 2, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p3_1d )
{
    // Compare static and dynamic P3 segment
    using convex_t = Simplex<1, 1, 1>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 3, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 3 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p1_2d )
{
    // Compare static and dynamic P1 triangle
    using convex_t = Simplex<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 1, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p2_2d )
{
    // Compare static and dynamic P2 triangle
    using convex_t = Simplex<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 2, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p3_2d )
{
    // Compare static and dynamic P3 triangle
    using convex_t = Simplex<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 3, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 3 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p4_2d )
{
    // Compare static and dynamic P4 triangle
    using convex_t = Simplex<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 4, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 4 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p1_3d )
{
    // Compare static and dynamic P1 tetrahedron
    using convex_t = Simplex<3, 1, 3>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 1, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p2_3d )
{
    // Compare static and dynamic P2 tetrahedron
    using convex_t = Simplex<3, 1, 3>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 2, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p3_3d )
{
    // Compare static and dynamic P3 tetrahedron
    using convex_t = Simplex<3, 1, 3>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 3, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 3 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_q1_2d )
{
    // Compare static and dynamic Q1 quadrilateral
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 1, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_q2_2d )
{
    // Compare static and dynamic Q2 quadrilateral
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 2, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_q3_2d )
{
    // Compare static and dynamic Q3 quadrilateral
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 3, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 3 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_q1_3d )
{
    // Compare static and dynamic Q1 hexahedron
    using convex_t = Hypercube<3, 1, 3>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 1, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 1 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_q2_3d )
{
    // Compare static and dynamic Q2 hexahedron
    using convex_t = Hypercube<3, 1, 3>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 2, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 2 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

//=============================================================================
// SECTION 4: Runtime Accessor Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_runtime_accessors_static )
{
    // Test runtime accessors for static order
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 3, double>;

    pointset_t pts;

    // Runtime accessors should return compile-time values for static order
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 3 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), pointset_t::numPoints );
    BOOST_CHECK_EQUAL( pts.runtimeNbPtsPerVertex(), pointset_t::nbPtsPerVertex );
    BOOST_CHECK_EQUAL( pts.runtimeNbPtsPerEdge(), pointset_t::nbPtsPerEdge );
    BOOST_CHECK_EQUAL( pts.runtimeNbPtsPerFace(), pointset_t::nbPtsPerFace );
}

BOOST_AUTO_TEST_CASE( test_pointset_runtime_accessors_dynamic )
{
    // Test runtime accessors for dynamic order
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    // P1
    pointset_t pts_p1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts_p1.runtimeOrder(), 1 );
    BOOST_CHECK_EQUAL( pts_p1.runtimeNumPoints(), 3 );
    BOOST_CHECK_EQUAL( pts_p1.runtimeNbPtsPerVertex(), 1 );
    BOOST_CHECK_EQUAL( pts_p1.runtimeNbPtsPerEdge(), 0 );
    BOOST_CHECK_EQUAL( pts_p1.runtimeNbPtsPerFace(), 0 );

    // P2
    pointset_t pts_p2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts_p2.runtimeOrder(), 2 );
    BOOST_CHECK_EQUAL( pts_p2.runtimeNumPoints(), 6 );
    BOOST_CHECK_EQUAL( pts_p2.runtimeNbPtsPerVertex(), 1 );
    BOOST_CHECK_EQUAL( pts_p2.runtimeNbPtsPerEdge(), 1 );
    BOOST_CHECK_EQUAL( pts_p2.runtimeNbPtsPerFace(), 0 );

    // P3
    pointset_t pts_p3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( pts_p3.runtimeOrder(), 3 );
    BOOST_CHECK_EQUAL( pts_p3.runtimeNumPoints(), 10 );
    BOOST_CHECK_EQUAL( pts_p3.runtimeNbPtsPerVertex(), 1 );
    BOOST_CHECK_EQUAL( pts_p3.runtimeNbPtsPerEdge(), 2 );
    BOOST_CHECK_EQUAL( pts_p3.runtimeNbPtsPerFace(), 1 );
}

//=============================================================================
// SECTION 5: Order 0 Tests (Special Case)
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_p0_static )
{
    // P0 elements have 1 centroid point
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 0, double>;

    static_assert( pointset_t::numPoints == 1, "P0 has 1 point" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 0 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 1 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 1 );
}

BOOST_AUTO_TEST_CASE( test_pointset_p0_dynamic )
{
    // Dynamic P0
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 0 ) );
    BOOST_CHECK_EQUAL( pts.runtimeOrder(), 0 );
    BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), 1 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 1 );
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_p0 )
{
    // Compare static and dynamic P0
    using convex_t = Simplex<2, 1, 2>;
    using pointset_static_t = PointSetEquiSpaced<convex_t, 0, double>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_static_t pts_static;
    pointset_dynamic_t pts_dynamic( RuntimeOrder( 0 ) );

    BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
    BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
    BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points() ) );
}

//=============================================================================
// SECTION 6: Multiple Dynamic Orders Test (Loop over orders)
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_multiple_orders_simplex_2d )
{
    // Test multiple orders in sequence for 2D simplex
    using convex_t = Simplex<2, 1, 2>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    for ( int order = 0; order <= 10; ++order )
    {
        // Use brace initialization to avoid most vexing parse
        pointset_dynamic_t pts{ RuntimeOrder{ static_cast<uint16_type>( order ) } };
        BOOST_CHECK_EQUAL( pts.runtimeOrder(), order );

        // Expected number of points: (order+1)*(order+2)/2
        uint32_type expected = ( order + 1 ) * ( order + 2 ) / 2;
        BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), expected );
        BOOST_CHECK_EQUAL( pts.points().size2(), expected );
    }
}

BOOST_AUTO_TEST_CASE( test_pointset_multiple_orders_hypercube_2d )
{
    // Test multiple orders in sequence for 2D hypercube
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    for ( int order = 0; order <= 10; ++order )
    {
        // Use brace initialization to avoid most vexing parse
        pointset_dynamic_t pts{ RuntimeOrder{ static_cast<uint16_type>( order ) } };
        BOOST_CHECK_EQUAL( pts.runtimeOrder(), order );

        // Expected number of points: (order+1)^2
        uint32_type expected = ( order + 1 ) * ( order + 1 );
        BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), expected );
        BOOST_CHECK_EQUAL( pts.points().size2(), expected );
    }
}

BOOST_AUTO_TEST_CASE( test_pointset_multiple_orders_simplex_3d )
{
    // Test multiple orders in sequence for 3D simplex
    using convex_t = Simplex<3, 1, 3>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    for ( int order = 0; order <= 10; ++order )
    {
        pointset_dynamic_t pts{ RuntimeOrder{ static_cast<uint16_type>( order ) } };
        BOOST_CHECK_EQUAL( pts.runtimeOrder(), order );

        // Expected number of points: (order+1)*(order+2)*(order+3)/6
        uint32_type expected = ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
        BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), expected );
        BOOST_CHECK_EQUAL( pts.points().size2(), expected );
    }
}

BOOST_AUTO_TEST_CASE( test_pointset_multiple_orders_hypercube_3d )
{
    // Test multiple orders in sequence for 3D hypercube
    using convex_t = Hypercube<3, 1, 3>;
    using pointset_dynamic_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    for ( int order = 0; order <= 10; ++order )
    {
        pointset_dynamic_t pts{ RuntimeOrder{ static_cast<uint16_type>( order ) } };
        BOOST_CHECK_EQUAL( pts.runtimeOrder(), order );

        // Expected number of points: (order+1)^3
        uint32_type expected = ( order + 1 ) * ( order + 1 ) * ( order + 1 );
        BOOST_CHECK_EQUAL( pts.runtimeNumPoints(), expected );
        BOOST_CHECK_EQUAL( pts.points().size2(), expected );
    }
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_simplex_2d_orders_0_to_10 )
{
    checkStaticDynamicPointsetRange<10, Simplex<2, 1, 2>>();
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_hypercube_2d_orders_0_to_10 )
{
    checkStaticDynamicPointsetRange<10, Hypercube<2, 1, 2>>();
}

BOOST_AUTO_TEST_SUITE_END()
