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
   \file test_dynamic_order_integration.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-05
   \brief Integration tests for type-level dynamic order
 */

#define BOOST_TEST_MODULE test_dynamic_order_integration
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/concepts.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelpoly/lagrange.hpp>

#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( dynamic_order_integration_suite )

//=============================================================================
// SECTION 1: Type System Foundation Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_dynamic_constant )
{
    // Dynamic constant should be -1 (Eigen-style)
    static_assert( Dynamic == -1, "Dynamic should equal -1" );
}

BOOST_AUTO_TEST_CASE( test_order_wrapper_basic )
{
    // RuntimeOrder wrapper for runtime order specification
    RuntimeOrder o1( 1 );
    RuntimeOrder o2( 2 );
    RuntimeOrder o3( 3 );

    BOOST_CHECK_EQUAL( o1.value, 1 );
    BOOST_CHECK_EQUAL( o2.value, 2 );
    BOOST_CHECK_EQUAL( o3.value, 3 );

    // Implicit conversion to uint16_type
    uint16_type val = o2;
    BOOST_CHECK_EQUAL( val, 2 );
}

BOOST_AUTO_TEST_CASE( test_type_distinctness_simplex )
{
    // Static and dynamic variants must be distinct types
    using simplex_p1 = Simplex<2, 1>;
    using simplex_p2 = Simplex<2, 2>;
    using simplex_dyn = Simplex<2, Dynamic>;

    static_assert( !std::is_same_v<simplex_p1, simplex_p2>,
                   "Different static orders should be different types" );
    static_assert( !std::is_same_v<simplex_p1, simplex_dyn>,
                   "Static and dynamic should be different types" );
    static_assert( !std::is_same_v<simplex_p2, simplex_dyn>,
                   "Static and dynamic should be different types" );
}

BOOST_AUTO_TEST_CASE( test_type_distinctness_lagrange )
{
    // Same for Lagrange
    using lagrange_p1 = Lagrange<1, Scalar>;
    using lagrange_p2 = Lagrange<2, Scalar>;
    using lagrange_dyn = Lagrange<Dynamic, Scalar>;

    static_assert( !std::is_same_v<lagrange_p1, lagrange_p2>,
                   "Different static orders should be different types" );
    static_assert( !std::is_same_v<lagrange_p1, lagrange_dyn>,
                   "Static and dynamic Lagrange should be different types" );
}

BOOST_AUTO_TEST_CASE( test_is_order_flags_simplex )
{
    // Static order types
    using simplex_static = Simplex<2, 1>;
    static_assert( simplex_static::is_order_static, "Expected static order" );
    static_assert( !simplex_static::is_order_dynamic, "Expected not dynamic" );

    // Dynamic order types
    using simplex_dynamic = Simplex<2, Dynamic>;
    static_assert( simplex_dynamic::is_order_dynamic, "Expected dynamic order" );
    static_assert( !simplex_dynamic::is_order_static, "Expected not static" );
}

BOOST_AUTO_TEST_CASE( test_is_order_flags_lagrange )
{
    // Same for Lagrange
    using lagrange_static = Lagrange<2, Scalar>;
    static_assert( lagrange_static::is_order_static, "Expected static order" );

    using lagrange_dynamic = Lagrange<Dynamic, Scalar>;
    static_assert( lagrange_dynamic::is_order_dynamic, "Expected dynamic order" );
}

BOOST_AUTO_TEST_CASE( test_static_only_accessors_guards )
{
    using simplex_static = Simplex<2, 2>;
    static_assert( requires { simplex_static::staticOrder(); } );
    static_assert( requires { simplex_static::staticNumPoints(); } );
    static_assert( simplex_static::staticOrder() == 2 );
    static_assert( simplex_static::staticNumPoints() == 6 );

    using hypercube_static = Hypercube<2, 2>;
    static_assert( requires { hypercube_static::staticOrder(); } );
    static_assert( requires { hypercube_static::staticNumPoints(); } );
    static_assert( hypercube_static::staticOrder() == 2 );
    static_assert( hypercube_static::staticNumPoints() == 9 );

    BOOST_CHECK( true );
}

//=============================================================================
// SECTION 2: Simplex Static vs Dynamic Comparison
//=============================================================================

BOOST_AUTO_TEST_CASE( test_simplex_static_vs_dynamic_values )
{
    // Compare values for P1, P2, P3 in 2D
    using simplex_2d_dyn = Simplex<2, Dynamic>;
    for ( uint16_type order = 1; order <= 10; ++order )
    {
        simplex_2d_dyn s_dyn{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( s_dyn.order(), order );

        // Compare with expected formula: (n+1)(n+2)/2
        uint16_type expected_points = ( order + 1 ) * ( order + 2 ) / 2;
        BOOST_CHECK_EQUAL( s_dyn.nPointsTotal(), expected_points );
    }
}

BOOST_AUTO_TEST_CASE( test_simplex_3d_static_vs_dynamic )
{
    // Compare 3D values
    using simplex_3d_dyn = Simplex<3, Dynamic>;
    for ( uint16_type order = 1; order <= 10; ++order )
    {
        simplex_3d_dyn s_dyn{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( s_dyn.order(), order );

        // Compare with expected formula: (n+1)(n+2)(n+3)/6
        uint16_type expected_points = ( order + 1 ) * ( order + 2 ) * ( order + 3 ) / 6;
        BOOST_CHECK_EQUAL( s_dyn.nPointsTotal(), expected_points );
    }
}

BOOST_AUTO_TEST_CASE( test_simplex_1d_static_vs_dynamic )
{
    // Compare 1D values
    using simplex_1d_dyn = Simplex<1, Dynamic>;
    for ( uint16_type order = 1; order <= 10; ++order )
    {
        simplex_1d_dyn s_dyn{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( s_dyn.order(), order );

        // 1D: n+1 points
        BOOST_CHECK_EQUAL( s_dyn.nPointsTotal(), order + 1 );
    }
}

//=============================================================================
// SECTION 3: Cross-validation with polyDims
//=============================================================================

BOOST_AUTO_TEST_CASE( test_polydims_matches_numpoints )
{
    // Verify polyDims matches numPoints for various orders
    using simplex_2d_dyn = Simplex<2, Dynamic>;
    using simplex_3d_dyn = Simplex<3, Dynamic>;
    using simplex_1d_dyn = Simplex<1, Dynamic>;
    using simplex_2d_p1 = Simplex<2, 1>;
    using simplex_3d_p1 = Simplex<3, 1>;
    using simplex_1d_p1 = Simplex<1, 1>;
    for ( uint16_type order = 0; order <= 10; ++order )
    {
        // 2D
        simplex_2d_dyn s2{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( s2.nPointsTotal(), simplex_2d_p1::polyDims( order ) );

        // 3D
        simplex_3d_dyn s3{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( s3.nPointsTotal(), simplex_3d_p1::polyDims( order ) );

        // 1D
        simplex_1d_dyn s1{ RuntimeOrder{ order } };
        BOOST_CHECK_EQUAL( s1.nPointsTotal(), simplex_1d_p1::polyDims( order ) );
    }
}

//=============================================================================
// SECTION 4: Per-entity Point Distribution Validation
//=============================================================================

BOOST_AUTO_TEST_CASE( test_point_distribution_2d )
{
    // Verify point distribution formula for 2D simplex
    // Total = 3*vertex + 3*edge + 1*face
    using simplex_2d_dyn = Simplex<2, Dynamic>;
    for ( uint16_type order = 1; order <= 10; ++order )
    {
        simplex_2d_dyn s{ RuntimeOrder{ order } };

        uint16_type pts_vertex = s.nPointsOnVertex();
        uint16_type pts_edge = s.nPointsOnEdge();
        uint16_type pts_face = s.nPointsOnFace();

        // 2D triangle: 3 vertices, 3 edges, 1 face (interior)
        uint16_type total = 3 * pts_vertex + 3 * pts_edge + pts_face;

        BOOST_CHECK_EQUAL( s.nPointsTotal(), total );
    }
}

BOOST_AUTO_TEST_CASE( test_point_distribution_3d )
{
    // Verify point distribution formula for 3D simplex
    // Total = 4*vertex + 6*edge + 4*face + 1*volume
    using simplex_3d_dyn = Simplex<3, Dynamic>;
    for ( uint16_type order = 1; order <= 10; ++order )
    {
        simplex_3d_dyn s{ RuntimeOrder{ order } };

        uint16_type pts_vertex = s.nPointsOnVertex();
        uint16_type pts_edge = s.nPointsOnEdge();
        uint16_type pts_face = s.nPointsOnFace();
        uint16_type pts_volume = s.nPointsOnVolume();

        // 3D tetra: 4 vertices, 6 edges, 4 faces, 1 volume
        uint16_type total = 4 * pts_vertex + 6 * pts_edge + 4 * pts_face + pts_volume;

        BOOST_CHECK_EQUAL( s.nPointsTotal(), total );
    }
}

//=============================================================================
// SECTION 5: Static Constants Match Dynamic Values
//=============================================================================

BOOST_AUTO_TEST_CASE( test_static_equals_dynamic_p1_2d )
{
    using static_t = Simplex<2, 1>;
    using dynamic_t = Simplex<2, Dynamic>;
    dynamic_t dyn{ RuntimeOrder{ 1 } };

    BOOST_CHECK_EQUAL( static_t::numPoints, dyn.nPointsTotal() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVertex, dyn.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerEdge, dyn.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerFace, dyn.nPointsOnFace() );
}

BOOST_AUTO_TEST_CASE( test_static_equals_dynamic_p2_2d )
{
    using static_t = Simplex<2, 2>;
    using dynamic_t = Simplex<2, Dynamic>;
    dynamic_t dyn{ RuntimeOrder{ 2 } };

    BOOST_CHECK_EQUAL( static_t::numPoints, dyn.nPointsTotal() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVertex, dyn.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerEdge, dyn.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerFace, dyn.nPointsOnFace() );
}

BOOST_AUTO_TEST_CASE( test_static_equals_dynamic_p3_2d )
{
    using static_t = Simplex<2, 3>;
    using dynamic_t = Simplex<2, Dynamic>;
    dynamic_t dyn{ RuntimeOrder{ 3 } };

    BOOST_CHECK_EQUAL( static_t::numPoints, dyn.nPointsTotal() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVertex, dyn.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerEdge, dyn.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerFace, dyn.nPointsOnFace() );
}

BOOST_AUTO_TEST_CASE( test_static_equals_dynamic_p1_3d )
{
    using static_t = Simplex<3, 1>;
    using dynamic_t = Simplex<3, Dynamic>;
    dynamic_t dyn{ RuntimeOrder{ 1 } };

    BOOST_CHECK_EQUAL( static_t::numPoints, dyn.nPointsTotal() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVertex, dyn.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerEdge, dyn.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerFace, dyn.nPointsOnFace() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVolume, dyn.nPointsOnVolume() );
}

BOOST_AUTO_TEST_CASE( test_static_equals_dynamic_p2_3d )
{
    using static_t = Simplex<3, 2>;
    using dynamic_t = Simplex<3, Dynamic>;
    dynamic_t dyn{ RuntimeOrder{ 2 } };

    BOOST_CHECK_EQUAL( static_t::numPoints, dyn.nPointsTotal() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVertex, dyn.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerEdge, dyn.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerFace, dyn.nPointsOnFace() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVolume, dyn.nPointsOnVolume() );
}

BOOST_AUTO_TEST_CASE( test_static_equals_dynamic_p10_2d )
{
    using static_t = Simplex<2, 10>;
    using dynamic_t = Simplex<2, Dynamic>;
    dynamic_t dyn{ RuntimeOrder{ 10 } };

    BOOST_CHECK_EQUAL( static_t::numPoints, dyn.nPointsTotal() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerVertex, dyn.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerEdge, dyn.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( static_t::nbPtsPerFace, dyn.nPointsOnFace() );
}

//=============================================================================
// SECTION 6: Lagrange Factory Integration
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_factory_static )
{
    // Verify Lagrange factory produces correct FE types for static orders
    using fe_p1 = typename Lagrange<1, Scalar>::template apply<2, 2>::type;
    using fe_p2 = typename Lagrange<2, Scalar>::template apply<2, 2>::type;
    using fe_p10 = typename Lagrange<10, Scalar>::template apply<2, 2>::type;

    fe_p1 p1;
    fe_p2 p2;
    fe_p10 p10;

    BOOST_CHECK_EQUAL( fe_p1::nOrder, 1 );
    BOOST_CHECK_EQUAL( fe_p2::nOrder, 2 );
    BOOST_CHECK_EQUAL( fe_p10::nOrder, 10 );
    BOOST_CHECK_EQUAL( fe_p1::nLocalDof, 3 );
    BOOST_CHECK_EQUAL( fe_p2::nLocalDof, 6 );
    BOOST_CHECK_EQUAL( fe_p10::nLocalDof, 66 );
}

BOOST_AUTO_TEST_CASE( test_lagrange_factory_dynamic )
{
    // Verify Lagrange<Dynamic> factory works
    using fe_dyn = typename Lagrange<Dynamic, Scalar>::template apply<2, 2>::type;

    fe_dyn fe;

    // Dynamic Lagrange uses placeholder order internally
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_SUITE_END()
