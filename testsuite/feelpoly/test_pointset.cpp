/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file test_pointset.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-06
   \brief Tests for PointSetEquiSpaced with static and dynamic order
 */

#define BOOST_TEST_MODULE test_pointset
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/equispaced.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feelpoly/gauss.hpp>
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

template<template<class, int, class> class PointSetType,
         int Order, typename ConvexType>
void checkStaticDynamicPointsetEquality(
    int interior = 0, double tol = 1e-12 )
{
    using pointset_static_t = PointSetType<ConvexType, Order, double>;
    using pointset_dynamic_t = PointSetType<ConvexType, Dynamic, double>;
    static_assert( HasPointCardinality<pointset_static_t> );
    static_assert( HasPointCardinality<pointset_dynamic_t> );

    pointset_static_t pts_static( interior );
    pointset_dynamic_t pts_dynamic{
        RuntimeOrder{ static_cast<uint16_type>( Order ) }, interior };

    BOOST_TEST_CONTEXT(
        "static/dynamic pointset equality order=" << Order <<
        " interior=" << interior )
    {
        BOOST_CHECK_EQUAL( pts_static.nPoints(), pts_dynamic.nPoints() );
        BOOST_CHECK_EQUAL( pts_static.points().size1(), pts_dynamic.points().size1() );
        BOOST_CHECK_EQUAL( pts_static.points().size2(), pts_dynamic.points().size2() );
        BOOST_CHECK( comparePointMatrices( pts_static.points(), pts_dynamic.points(), tol ) );

        using reference_type = Reference<
            ConvexType, ConvexType::nDim, ConvexType::nOrder,
            ConvexType::nRealDim, double>;
        reference_type reference;
        for ( uint16_type entityDim = 0;
              entityDim <= ConvexType::nDim; ++entityDim )
        {
            for ( int entity = reference.entityRange( entityDim ).begin();
                  entity < reference.entityRange( entityDim ).end();
                  ++entity )
            {
                auto staticRestriction = pts_static.pointsBySubEntity(
                    entityDim, entity, false );
                auto dynamicRestriction = pts_dynamic.pointsBySubEntity(
                    entityDim, entity, false );
                BOOST_CHECK_EQUAL(
                    staticRestriction.size2(),
                    dynamicRestriction.size2() );
                BOOST_CHECK( comparePointMatrices(
                    staticRestriction, dynamicRestriction, tol ) );
            }
        }
    }
}

template<template<class, int, class> class PointSetType,
         int MaxOrder, typename ConvexType, int CurrentOrder = 0>
void checkStaticDynamicPointsetRange(
    int interior = 0, double tol = 1e-12 )
{
    checkStaticDynamicPointsetEquality<
        PointSetType, CurrentOrder, ConvexType>( interior, tol );
    if constexpr ( CurrentOrder < MaxOrder )
        checkStaticDynamicPointsetRange<
            PointSetType, MaxOrder, ConvexType, CurrentOrder+1>(
                interior, tol );
}

template<template<class, int, class> class PointSetType,
         typename ConvexType, typename Cardinality>
void checkStoredPointCountRange()
{
    using pointset_type =
        PointSetType<ConvexType, Dynamic, double>;
    static_assert( HasPointCardinality<pointset_type> );

    for ( uint16_type order = 0; order <= 10; ++order )
    {
        pointset_type complete{ RuntimeOrder{ order }, 0 };
        pointset_type interior{ RuntimeOrder{ order }, 1 };
        const auto fullCount = Cardinality::total(
            ConvexType::nDim, order );
        const auto interiorCount = Cardinality::onEntity(
            ConvexType::nDim, ConvexType::nDim, order );

        BOOST_TEST_CONTEXT( "stored/full cardinality order=" << order )
        {
            BOOST_CHECK_EQUAL( complete.order(), order );
            BOOST_CHECK_EQUAL( interior.order(), order );
            BOOST_CHECK_EQUAL( complete.nPoints(), fullCount );
            BOOST_CHECK_EQUAL( complete.nPoints(), complete.points().size2() );
            BOOST_CHECK_EQUAL( interior.nPoints(), interiorCount );
            BOOST_CHECK_EQUAL( interior.nPoints(), interior.points().size2() );

            auto completeInterior = complete.pointsBySubEntity(
                ConvexType::nDim, 0, false );
            auto interiorRestriction = interior.pointsBySubEntity(
                ConvexType::nDim, 0, false );
            BOOST_CHECK_EQUAL( completeInterior.size2(), interiorCount );
            BOOST_CHECK_EQUAL( interiorRestriction.size2(), interiorCount );
            BOOST_CHECK( comparePointMatrices(
                completeInterior, interiorRestriction ) );

            if constexpr ( requires {
                               complete.nPointsOnVertex();
                               complete.nPointsOnEdge();
                               complete.nPointsOnFace();
                               complete.nPointsOnVolume();
                           } )
            {
                BOOST_CHECK_EQUAL(
                    complete.nPointsOnVertex(),
                    Cardinality::perVertex( ConvexType::nDim, order ) );
                BOOST_CHECK_EQUAL(
                    complete.nPointsOnEdge(),
                    Cardinality::perEdge( ConvexType::nDim, order ) );
                BOOST_CHECK_EQUAL(
                    complete.nPointsOnFace(),
                    Cardinality::perFace( ConvexType::nDim, order ) );
                BOOST_CHECK_EQUAL(
                    complete.nPointsOnVolume(),
                    Cardinality::perVolume( ConvexType::nDim, order ) );
            }
        }
    }
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
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 2 );
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
    BOOST_CHECK_EQUAL( pts.order(), 2 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p1_2d )
{
    // Test P1 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 1, double>;

    static_assert( pointset_t::numPoints == 3, "P1 triangle has 3 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 3 );
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
    BOOST_CHECK_EQUAL( pts.order(), 2 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 6 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 6 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p3_2d )
{
    // Test P3 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 3, double>;

    static_assert( pointset_t::numPoints == 10, "P3 triangle has 10 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.order(), 3 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 10 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_p1_3d )
{
    // Test P1 tetrahedron (3D)
    using convex_t = Simplex<3, 1, 3>;
    using pointset_t = PointSetEquiSpaced<convex_t, 1, double>;

    static_assert( pointset_t::numPoints == 4, "P1 tetra has 4 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 4 );
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
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 4 );
}

BOOST_AUTO_TEST_CASE( test_pointset_static_q2_2d )
{
    // Test Q2 quadrilateral (2D hypercube)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 2, double>;

    static_assert( pointset_t::numPoints == 9, "Q2 quad has 9 points" );

    pointset_t pts;
    BOOST_CHECK_EQUAL( pts.order(), 2 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 9 );
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
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 2 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 2 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p2_1d )
{
    // Test dynamic P2 segment (1D)
    using convex_t = Simplex<1, 1, 1>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts.order(), 2 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p1_2d )
{
    // Test dynamic P1 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 3 );
    BOOST_CHECK_EQUAL( pts.points().size1(), 2 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 3 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p2_2d )
{
    // Test dynamic P2 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts.order(), 2 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 6 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p3_2d )
{
    // Test dynamic P3 triangle (2D)
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( pts.order(), 3 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 10 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_p1_3d )
{
    // Test dynamic P1 tetrahedron (3D)
    using convex_t = Simplex<3, 1, 3>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 4 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_q1_2d )
{
    // Test dynamic Q1 quadrilateral (2D)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts.order(), 1 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 4 );
}

BOOST_AUTO_TEST_CASE( test_pointset_dynamic_q2_2d )
{
    // Test dynamic Q2 quadrilateral (2D)
    using convex_t = Hypercube<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts.order(), 2 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 9 );
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
// SECTION 4: Semantic Cardinality Accessor Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_pointset_semantic_accessors_static )
{
    // Static and dynamic point sets expose the same object-level vocabulary.
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, 3, double>;

    pointset_t pts;

    // Static objects report their effective values through the same API.
    BOOST_CHECK_EQUAL( pts.order(), 3 );
    BOOST_CHECK_EQUAL( pts.nPoints(), pointset_t::numPoints );
    BOOST_CHECK_EQUAL( pts.nPointsOnVertex(), pointset_t::nbPtsPerVertex );
    BOOST_CHECK_EQUAL( pts.nPointsOnEdge(), pointset_t::nbPtsPerEdge );
    BOOST_CHECK_EQUAL( pts.nPointsOnFace(), pointset_t::nbPtsPerFace );
}

BOOST_AUTO_TEST_CASE( test_pointset_semantic_accessors_dynamic )
{
    // Dynamic objects report values computed from their stored effective order.
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    // P1
    pointset_t pts_p1( RuntimeOrder( 1 ) );
    BOOST_CHECK_EQUAL( pts_p1.order(), 1 );
    BOOST_CHECK_EQUAL( pts_p1.nPoints(), 3 );
    BOOST_CHECK_EQUAL( pts_p1.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( pts_p1.nPointsOnEdge(), 0 );
    BOOST_CHECK_EQUAL( pts_p1.nPointsOnFace(), 0 );

    // P2
    pointset_t pts_p2( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( pts_p2.order(), 2 );
    BOOST_CHECK_EQUAL( pts_p2.nPoints(), 6 );
    BOOST_CHECK_EQUAL( pts_p2.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( pts_p2.nPointsOnEdge(), 1 );
    BOOST_CHECK_EQUAL( pts_p2.nPointsOnFace(), 0 );

    // P3
    pointset_t pts_p3( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( pts_p3.order(), 3 );
    BOOST_CHECK_EQUAL( pts_p3.nPoints(), 10 );
    BOOST_CHECK_EQUAL( pts_p3.nPointsOnVertex(), 1 );
    BOOST_CHECK_EQUAL( pts_p3.nPointsOnEdge(), 2 );
    BOOST_CHECK_EQUAL( pts_p3.nPointsOnFace(), 1 );
}

BOOST_AUTO_TEST_CASE( test_pointset_quadrature_interpolation_cardinality_contract )
{
    using interpolation_type =
        PointSetEquiSpaced<Simplex<2, 1, 2>, 3, double>;
    using quadrature_type = Gauss<Simplex<2, 1>, 4, double>;

    static_assert( HasPointCardinality<interpolation_type> );
    static_assert( HasPointCardinality<quadrature_type> );

    interpolation_type interpolation;
    quadrature_type quadrature;

    BOOST_CHECK_EQUAL(
        interpolation.nPoints(), interpolation.points().size2() );
    BOOST_CHECK_EQUAL(
        quadrature.nPoints(), quadrature.points().size2() );
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
    BOOST_CHECK_EQUAL( pts.order(), 0 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 1 );
    BOOST_CHECK_EQUAL( pts.points().size2(), 1 );
}

BOOST_AUTO_TEST_CASE( test_pointset_p0_dynamic )
{
    // Dynamic P0
    using convex_t = Simplex<2, 1, 2>;
    using pointset_t = PointSetEquiSpaced<convex_t, Dynamic, double>;

    pointset_t pts( RuntimeOrder( 0 ) );
    BOOST_CHECK_EQUAL( pts.order(), 0 );
    BOOST_CHECK_EQUAL( pts.nPoints(), 1 );
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
        BOOST_CHECK_EQUAL( pts.order(), order );

        // Expected number of points: (order+1)*(order+2)/2
        uint32_type expected = ( order + 1 ) * ( order + 2 ) / 2;
        BOOST_CHECK_EQUAL( pts.nPoints(), expected );
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
        BOOST_CHECK_EQUAL( pts.order(), order );

        // Expected number of points: (order+1)^2
        uint32_type expected = ( order + 1 ) * ( order + 1 );
        BOOST_CHECK_EQUAL( pts.nPoints(), expected );
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
        BOOST_CHECK_EQUAL( pts.order(), order );

        // Expected number of points: (order+1)*(order+2)*(order+3)/6
        uint32_type expected = ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
        BOOST_CHECK_EQUAL( pts.nPoints(), expected );
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
        BOOST_CHECK_EQUAL( pts.order(), order );

        // Expected number of points: (order+1)^3
        uint32_type expected = ( order + 1 ) * ( order + 1 ) * ( order + 1 );
        BOOST_CHECK_EQUAL( pts.nPoints(), expected );
        BOOST_CHECK_EQUAL( pts.points().size2(), expected );
    }
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_simplex_2d_orders_0_to_10 )
{
    checkStaticDynamicPointsetRange<
        PointSetEquiSpaced, 10, Simplex<2, 1, 2>>();
}

BOOST_AUTO_TEST_CASE( test_pointset_compare_hypercube_2d_orders_0_to_10 )
{
    checkStaticDynamicPointsetRange<
        PointSetEquiSpaced, 10, Hypercube<2, 1, 2>>();
}

BOOST_AUTO_TEST_CASE( test_pointset_interior_static_dynamic_orders_0_to_10 )
{
    checkStaticDynamicPointsetRange<
        PointSetEquiSpaced, 10, Simplex<2, 1, 2>>( 1 );
    checkStaticDynamicPointsetRange<
        PointSetEquiSpaced, 10, Hypercube<2, 1, 2>>( 1 );
}

BOOST_AUTO_TEST_CASE( test_pointset_stored_and_full_cardinality_orders_0_to_10 )
{
    checkStoredPointCountRange<
        PointSetEquiSpaced, Simplex<2, 1, 2>,
        Feel::detail::SimplexCardinality>();
    checkStoredPointCountRange<
        PointSetEquiSpaced, Hypercube<2, 1, 2>,
        Feel::detail::HypercubeCardinality>();
    checkStoredPointCountRange<
        PointSetEquiSpaced, Simplex<3, 1, 3>,
        Feel::detail::SimplexCardinality>();
    checkStoredPointCountRange<
        PointSetEquiSpaced, Hypercube<3, 1, 3>,
        Feel::detail::HypercubeCardinality>();

    checkStoredPointCountRange<
        PointSetFekete, Simplex<2, 1, 2>,
        Feel::detail::SimplexCardinality>();
    checkStoredPointCountRange<
        PointSetFekete, Hypercube<2, 1, 2>,
        Feel::detail::HypercubeCardinality>();
}

BOOST_AUTO_TEST_CASE( test_fekete_static_dynamic_orders_0_to_10 )
{
    checkStaticDynamicPointsetRange<
        PointSetFekete, 10, Simplex<2, 1, 2>>();
    checkStaticDynamicPointsetRange<
        PointSetFekete, 10, Hypercube<2, 1, 2>>();
    checkStaticDynamicPointsetRange<
        PointSetFekete, 10, Simplex<2, 1, 2>>( 1 );
    checkStaticDynamicPointsetRange<
        PointSetFekete, 10, Hypercube<2, 1, 2>>( 1 );
}

BOOST_AUTO_TEST_SUITE_END()
