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
   \file test_lagrange_dynamic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-05
   \brief Tests for Lagrange<Dynamic> with type-level dynamic polynomial order
 */

#define BOOST_TEST_MODULE test_lagrange_dynamic
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/simplex.hpp>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <string>
#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename FEType>
void checkKroneckerDeltaProperty( FEType const& fe, double tol, std::string const& label )
{
    auto const pts = fe.points();
    auto const eval = fe.evaluate( pts );

    BOOST_REQUIRE_EQUAL( eval.size1(), eval.size2() );

    double maxDiagError = 0.0;
    double maxOffDiag = 0.0;
    for ( std::size_t i = 0; i < eval.size1(); ++i )
    {
        for ( std::size_t j = 0; j < eval.size2(); ++j )
        {
            const double v = static_cast<double>( eval( i, j ) );
            if ( i == j )
                maxDiagError = std::max( maxDiagError, std::abs( v - 1.0 ) );
            else
                maxOffDiag = std::max( maxOffDiag, std::abs( v ) );
        }
    }

    BOOST_TEST_CONTEXT( "Kronecker delta check: " << label )
    {
        BOOST_CHECK_SMALL( maxDiagError, tol );
        BOOST_CHECK_SMALL( maxOffDiag, tol );
    }
}

template <typename FEType>
void checkPartitionOfUnity( FEType const& fe,
                            ublas::matrix<double> const& pts,
                            double tol,
                            std::string const& label )
{
    auto const eval = fe.evaluate( pts );
    double maxError = 0.0;
    for ( std::size_t j = 0; j < eval.size2(); ++j )
    {
        double sumPhi = 0.0;
        for ( std::size_t i = 0; i < eval.size1(); ++i )
            sumPhi += static_cast<double>( eval( i, j ) );
        maxError = std::max( maxError, std::abs( sumPhi - 1.0 ) );
    }

    BOOST_TEST_CONTEXT( "Partition of unity: " << label )
    {
        BOOST_CHECK_SMALL( maxError, tol );
    }
}

template <typename MatrixA, typename MatrixB>
void checkMatricesNear( MatrixA const& a,
                        MatrixB const& b,
                        double tol,
                        std::string const& label )
{
    BOOST_TEST_CONTEXT( label )
    {
        BOOST_REQUIRE_EQUAL( a.size1(), b.size1() );
        BOOST_REQUIRE_EQUAL( a.size2(), b.size2() );

        double maxAbsDiff = 0.0;
        for ( std::size_t i = 0; i < a.size1(); ++i )
        {
            for ( std::size_t j = 0; j < a.size2(); ++j )
            {
                const double diff = std::abs( static_cast<double>( a( i, j ) ) -
                                              static_cast<double>( b( i, j ) ) );
                maxAbsDiff = std::max( maxAbsDiff, diff );
            }
        }
        BOOST_CHECK_SMALL( maxAbsDiff, tol );
    }
}

template <int Dim, int Order>
void checkStaticRuntimeLagrangeAgreement( double tol )
{
    using lagrange_static_t = Lagrange<Order, Scalar, Continuous, PointSetEquiSpaced>;
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar, Continuous, PointSetEquiSpaced>;
    using fe_static_t = typename lagrange_static_t::template apply<Dim, Dim, double, Simplex<Dim>>::type;
    using fe_dynamic_t = typename lagrange_dynamic_t::template apply<Dim, Dim, double, Simplex<Dim>>::type;

    fe_static_t fe_static;
    fe_dynamic_t fe_dynamic{ RuntimeOrder{ Order } };

    auto const pts_static = fe_static.points();
    auto const pts_dynamic = fe_dynamic.points();

    const std::string prefix = "Dim=" + std::to_string( Dim ) + " P" + std::to_string( Order ) + " ";

    checkMatricesNear( pts_static, pts_dynamic, tol, prefix + "node coordinates" );
    checkMatricesNear( fe_static.evaluate( pts_static ),
                       fe_dynamic.evaluate( pts_static ),
                       tol,
                       prefix + "basis eval on static nodes" );
    checkMatricesNear( fe_static.evaluate( pts_dynamic ),
                       fe_dynamic.evaluate( pts_dynamic ),
                       tol,
                       prefix + "basis eval on dynamic nodes" );

    ublas::matrix<double> probePts( Dim, 4 );
    const std::size_t nPts = pts_static.size2();
    BOOST_REQUIRE_GT( nPts, std::size_t( 0 ) );
    for ( std::size_t j = 0; j < probePts.size2(); ++j )
    {
        const std::size_t j0 = j % nPts;
        const std::size_t j1 = ( j + 1 ) % nPts;
        const std::size_t j2 = ( j + 2 ) % nPts;
        for ( std::size_t d = 0; d < probePts.size1(); ++d )
            probePts( d, j ) = 0.5 * pts_static( d, j0 ) + 0.3 * pts_static( d, j1 ) + 0.2 * pts_static( d, j2 );
    }

    checkMatricesNear( fe_static.evaluate( probePts ),
                       fe_dynamic.evaluate( probePts ),
                       tol,
                       prefix + "basis eval on probe points" );
}

template <int Dim, int Order, typename ConvexType>
void checkStaticRuntimeLagrangeAgreementDefaultPointSet( double tol, std::string const& labelPrefix )
{
    using lagrange_static_t = Lagrange<Order, Scalar, Continuous>;
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar, Continuous>;
    using fe_static_t = typename lagrange_static_t::template apply<Dim, Dim, double, ConvexType>::type;
    using fe_dynamic_t = typename lagrange_dynamic_t::template apply<Dim, Dim, double, ConvexType>::type;

    fe_static_t fe_static;
    fe_dynamic_t fe_dynamic{ RuntimeOrder{ Order } };

    auto const pts_static = fe_static.points();
    auto const pts_dynamic = fe_dynamic.points();

    const std::string prefix = labelPrefix + " Dim=" + std::to_string( Dim ) + " P" + std::to_string( Order ) + " ";

    checkMatricesNear( pts_static, pts_dynamic, tol, prefix + "node coordinates" );
    checkMatricesNear( fe_static.evaluate( pts_static ),
                       fe_dynamic.evaluate( pts_static ),
                       tol,
                       prefix + "basis eval on static nodes" );
    checkMatricesNear( fe_static.evaluate( pts_dynamic ),
                       fe_dynamic.evaluate( pts_dynamic ),
                       tol,
                       prefix + "basis eval on dynamic nodes" );
}
}

BOOST_AUTO_TEST_SUITE( lagrange_dynamic_suite )

//=============================================================================
// SECTION 1: Backward Compatibility Tests
// Ensure existing static-order Lagrange still works exactly as before
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_static_order_backward_compat )
{
    // Test that Lagrange<1, Scalar> (P1 scalar) works as before
    using lagrange_p1_t = Lagrange<1, Scalar>;

    static_assert( lagrange_p1_t::nOrder == 1, "Expected order 1" );
    static_assert( lagrange_p1_t::is_order_static, "Expected static order" );
    static_assert( !lagrange_p1_t::is_order_dynamic, "Expected not dynamic" );

    // Apply to get the actual finite element type
    using fe_t = typename lagrange_p1_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_t::nDim == 2, "Expected 2D" );
    static_assert( fe_t::nOrder == 1, "Expected order 1" );
    static_assert( fe_t::nLocalDof == 3, "P1 triangle has 3 DOFs" );
    static_assert( fe_t::is_scalar, "Expected scalar field" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_order_p2 )
{
    // Test P2 Lagrange
    using lagrange_p2_t = Lagrange<2, Scalar>;

    static_assert( lagrange_p2_t::nOrder == 2, "Expected order 2" );
    static_assert( lagrange_p2_t::is_order_static, "Expected static order" );

    using fe_t = typename lagrange_p2_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_t::nOrder == 2, "Expected order 2" );
    static_assert( fe_t::nLocalDof == 6, "P2 triangle has 6 DOFs" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_order_p10 )
{
    // Test P10 static Lagrange as high-order compile-time coverage
    using lagrange_p10_t = Lagrange<10, Scalar>;
    using fe_t = typename lagrange_p10_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_t::nOrder == 10, "Expected order 10" );
    static_assert( fe_t::nLocalDof == 66, "P10 triangle has 66 DOFs" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_order_vectorial )
{
    // Test vectorial Lagrange (P1 vector field)
    using lagrange_vec_t = Lagrange<1, Vectorial>;

    using fe_t = typename lagrange_vec_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_t::is_vectorial, "Expected vectorial field" );
    static_assert( fe_t::nComponents == 2, "2D vector has 2 components" );
    static_assert( fe_t::nLocalDof == 3, "3 DOFs per component" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_order_3d )
{
    // Test 3D Lagrange (P1 on tetrahedron)
    using lagrange_p1_t = Lagrange<1, Scalar>;

    using fe_t = typename lagrange_p1_t::template apply<3, 3, double, Simplex<3>>::type;

    static_assert( fe_t::nDim == 3, "Expected 3D" );
    static_assert( fe_t::nLocalDof == 4, "P1 tetra has 4 DOFs" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_kronecker_delta_property )
{
    using fe_p1_t = typename Lagrange<1, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using fe_p2_t = typename Lagrange<2, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using fe_p3_t = typename Lagrange<3, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    using fe_p10_t = typename Lagrange<10, Scalar>::template apply<2, 2, double, Simplex<2>>::type;

    checkKroneckerDeltaProperty( fe_p1_t{}, 1e-12, "static P1" );
    checkKroneckerDeltaProperty( fe_p2_t{}, 1e-11, "static P2" );
    checkKroneckerDeltaProperty( fe_p3_t{}, 1e-10, "static P3" );
    checkKroneckerDeltaProperty( fe_p10_t{}, 1e-8, "static P10" );
}

//=============================================================================
// SECTION 2: Type-Level Dynamic Order Tests
// Test Lagrange<Dynamic> type traits
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_type_traits )
{
    // Test that Lagrange<Dynamic> is a distinct type with correct traits
    using lagrange_static_t = Lagrange<2, Scalar>;
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar>;

    // Types should be different
    static_assert( !std::is_same_v<lagrange_static_t, lagrange_dynamic_t>,
                   "Static and dynamic Lagrange must be different types" );

    // Dynamic Lagrange should have is_order_dynamic = true
    static_assert( lagrange_dynamic_t::is_order_dynamic,
                   "Expected is_order_dynamic true for Lagrange<Dynamic>" );
    static_assert( !lagrange_dynamic_t::is_order_static,
                   "Expected is_order_static false for Lagrange<Dynamic>" );

    // Static Lagrange should have is_order_static = true
    static_assert( lagrange_static_t::is_order_static,
                   "Expected is_order_static true for Lagrange<2>" );
    static_assert( !lagrange_static_t::is_order_dynamic,
                   "Expected is_order_dynamic false for Lagrange<2>" );

    // nOrder should be 0 for dynamic (runtime order not known at compile time)
    static_assert( lagrange_dynamic_t::nOrder == 0,
                   "Expected nOrder == 0 for dynamic" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_apply )
{
    // Test that apply metafunction works for dynamic Lagrange
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar>;

    // Get the FE type - for dynamic, this uses a placeholder order internally
    using fe_type = typename lagrange_dynamic_t::template apply<2, 2, double, Simplex<2>>::type;

    // The FE type exists and compiles
    static_assert( fe_type::nDim == 2, "Expected 2D" );
    static_assert( fe_type::is_scalar, "Expected scalar field" );

    // We can construct it (uses placeholder order)
    fe_type fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_vectorial )
{
    // Test vectorial dynamic Lagrange
    using lagrange_dynamic_vec_t = Lagrange<Dynamic, Vectorial>;

    static_assert( lagrange_dynamic_vec_t::is_order_dynamic,
                   "Expected dynamic order" );

    using fe_t = typename lagrange_dynamic_vec_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_t::is_vectorial, "Expected vectorial field" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_3d )
{
    // Test 3D dynamic Lagrange
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar>;

    using fe_t = typename lagrange_dynamic_t::template apply<3, 3, double, Simplex<3>>::type;

    static_assert( fe_t::nDim == 3, "Expected 3D" );

    fe_t fe;
    BOOST_CHECK_EQUAL( fe.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_runtime_orders_up_to_10 )
{
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar>;
    using fe_dyn_t = typename lagrange_dynamic_t::template apply<2, 2, double, Simplex<2>>::type;
    using simplex_2d_p1 = Simplex<2, 1>;

    for ( uint16_type order = 1; order <= 10; ++order )
    {
        fe_dyn_t fe{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( fe.runtimeOrder(), order );
        BOOST_CHECK_EQUAL( fe.runtimeLocalDof(), simplex_2d_p1::polyDims( order ) );

        const auto dofPerVertex = fe.runtimeDofPerVertex();
        const auto dofPerEdge = fe.runtimeDofPerEdge();
        const auto dofPerFace = fe.runtimeDofPerFace();
        BOOST_CHECK_EQUAL( fe.runtimeLocalDof(), 3 * dofPerVertex + 3 * dofPerEdge + dofPerFace );
    }
}

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_kronecker_delta_property )
{
    using lagrange_dynamic_t = Lagrange<Dynamic, Scalar>;
    using fe_dyn_t = typename lagrange_dynamic_t::template apply<2, 2, double, Simplex<2>>::type;

    for ( uint16_type order : { uint16_type( 1 ), uint16_type( 2 ), uint16_type( 3 ), uint16_type( 10 ) } )
    {
        fe_dyn_t fe{ RuntimeOrder{ order } };
        const double tol = ( order >= 10 ) ? 1e-8 : 1e-10;
        checkKroneckerDeltaProperty( fe, tol, "dynamic P" + std::to_string( order ) );
    }
}

BOOST_AUTO_TEST_CASE( test_lagrange_dynamic_equispaced_policy )
{
    // Explicitly request PointSetEquiSpaced policy and verify runtime-order behavior.
    using lagrange_dynamic_equi_t = Lagrange<Dynamic, Scalar, Continuous, PointSetEquiSpaced>;
    using fe_dyn_equi_t = typename lagrange_dynamic_equi_t::template apply<2, 2, double, Simplex<2>>::type;
    using simplex_2d_p1 = Simplex<2, 1>;

    static_assert( lagrange_dynamic_equi_t::is_order_dynamic, "Expected dynamic order" );

    for ( uint16_type order : { uint16_type( 1 ), uint16_type( 2 ), uint16_type( 3 ), uint16_type( 10 ) } )
    {
        fe_dyn_equi_t fe{ RuntimeOrder{ order } };

        BOOST_CHECK_EQUAL( fe.runtimeOrder(), order );
        BOOST_CHECK_EQUAL( fe.runtimeLocalDof(), simplex_2d_p1::polyDims( order ) );

        const double tol = ( order >= 10 ) ? 1e-8 : 1e-10;
        checkKroneckerDeltaProperty( fe, tol, "dynamic equispaced P" + std::to_string( order ) );
    }
}

BOOST_AUTO_TEST_CASE( test_lagrange_partition_of_unity )
{
    ublas::matrix<double> pts( 2, 4 );
    pts( 0, 0 ) = -0.8; pts( 1, 0 ) = -0.8;
    pts( 0, 1 ) =  0.0; pts( 1, 1 ) = -0.8;
    pts( 0, 2 ) = -0.2; pts( 1, 2 ) = -0.1;
    pts( 0, 3 ) = -0.6; pts( 1, 3 ) =  0.2;

    checkPartitionOfUnity(
        typename Lagrange<1, Scalar>::template apply<2, 2, double, Simplex<2>>::type{},
        pts, 1e-12, "static P1" );
    checkPartitionOfUnity(
        typename Lagrange<2, Scalar>::template apply<2, 2, double, Simplex<2>>::type{},
        pts, 1e-11, "static P2" );
    checkPartitionOfUnity(
        typename Lagrange<3, Scalar>::template apply<2, 2, double, Simplex<2>>::type{},
        pts, 1e-10, "static P3" );
    checkPartitionOfUnity(
        typename Lagrange<10, Scalar>::template apply<2, 2, double, Simplex<2>>::type{},
        pts, 1e-8, "static P10" );

    using fe_dyn_t = typename Lagrange<Dynamic, Scalar>::template apply<2, 2, double, Simplex<2>>::type;
    checkPartitionOfUnity( fe_dyn_t{ RuntimeOrder{ 1 } }, pts, 1e-12, "dynamic P1" );
    checkPartitionOfUnity( fe_dyn_t{ RuntimeOrder{ 2 } }, pts, 1e-11, "dynamic P2" );
    checkPartitionOfUnity( fe_dyn_t{ RuntimeOrder{ 3 } }, pts, 1e-10, "dynamic P3" );
    checkPartitionOfUnity( fe_dyn_t{ RuntimeOrder{ 10 } }, pts, 1e-8, "dynamic P10" );
}

//=============================================================================
// SECTION 3: Comparison Tests - Static vs Dynamic Types
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_type_distinctness )
{
    // Verify that different orders produce different types
    using lagrange_p1 = Lagrange<1, Scalar>;
    using lagrange_p2 = Lagrange<2, Scalar>;
    using lagrange_p3 = Lagrange<3, Scalar>;
    using lagrange_dyn = Lagrange<Dynamic, Scalar>;

    static_assert( !std::is_same_v<lagrange_p1, lagrange_p2>,
                   "P1 and P2 should be different types" );
    static_assert( !std::is_same_v<lagrange_p1, lagrange_p3>,
                   "P1 and P3 should be different types" );
    static_assert( !std::is_same_v<lagrange_p1, lagrange_dyn>,
                   "P1 and Dynamic should be different types" );
    static_assert( !std::is_same_v<lagrange_p2, lagrange_dyn>,
                   "P2 and Dynamic should be different types" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_fe_static_values )
{
    // Verify that static FE values are correct for different orders
    using fe_p1 = typename Lagrange<1, Scalar>::template apply<2, 2>::type;
    using fe_p2 = typename Lagrange<2, Scalar>::template apply<2, 2>::type;
    using fe_p3 = typename Lagrange<3, Scalar>::template apply<2, 2>::type;

    // P1 triangle: 3 vertices = 3 DOFs
    static_assert( fe_p1::nLocalDof == 3, "P1 has 3 DOFs" );
    static_assert( fe_p1::nOrder == 1, "P1 has order 1" );

    // P2 triangle: 3 vertices + 3 edge midpoints = 6 DOFs
    static_assert( fe_p2::nLocalDof == 6, "P2 has 6 DOFs" );
    static_assert( fe_p2::nOrder == 2, "P2 has order 2" );

    // P3 triangle: 3 + 6 + 1 = 10 DOFs
    static_assert( fe_p3::nLocalDof == 10, "P3 has 10 DOFs" );
    static_assert( fe_p3::nOrder == 3, "P3 has order 3" );

    // Verify at runtime too
    fe_p1 p1;
    fe_p2 p2;
    fe_p3 p3;

    BOOST_CHECK_EQUAL( fe_p1::nLocalDof, 3 );
    BOOST_CHECK_EQUAL( fe_p2::nLocalDof, 6 );
    BOOST_CHECK_EQUAL( fe_p3::nLocalDof, 10 );
}

//=============================================================================
// SECTION 4: Continuity Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_continuity_static )
{
    // Test continuous vs discontinuous static Lagrange
    using lagrange_cont_t = Lagrange<2, Scalar, Continuous>;
    using lagrange_disc_t = Lagrange<2, Scalar, Discontinuous>;

    using fe_cont_t = typename lagrange_cont_t::template apply<2, 2, double, Simplex<2>>::type;
    using fe_disc_t = typename lagrange_disc_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_cont_t::isContinuous, "Expected continuous" );
    static_assert( !fe_disc_t::isContinuous, "Expected discontinuous" );

    fe_cont_t fe_cont;
    fe_disc_t fe_disc;

    BOOST_CHECK_EQUAL( fe_cont.familyName(), "lagrange" );
    BOOST_CHECK_EQUAL( fe_disc.familyName(), "lagrange" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_continuity_dynamic )
{
    // Test continuous vs discontinuous dynamic Lagrange
    using lagrange_cont_t = Lagrange<Dynamic, Scalar, Continuous>;
    using lagrange_disc_t = Lagrange<Dynamic, Scalar, Discontinuous>;

    static_assert( lagrange_cont_t::is_order_dynamic, "Expected dynamic" );
    static_assert( lagrange_disc_t::is_order_dynamic, "Expected dynamic" );

    using fe_cont_t = typename lagrange_cont_t::template apply<2, 2, double, Simplex<2>>::type;
    using fe_disc_t = typename lagrange_disc_t::template apply<2, 2, double, Simplex<2>>::type;

    static_assert( fe_cont_t::isContinuous, "Expected continuous" );
    static_assert( !fe_disc_t::isContinuous, "Expected discontinuous" );
}

//=============================================================================
// SECTION 5: ChangeTag Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_change_tag_static )
{
    using lagrange_t = Lagrange<2, Scalar>;
    using lagrange_tagged_t = typename lagrange_t::template ChangeTag<1>::type;

    static_assert( lagrange_t::TAG == 0, "Original has TAG 0" );
    static_assert( lagrange_tagged_t::TAG == 1, "Tagged has TAG 1" );
    static_assert( lagrange_tagged_t::nOrder == 2, "Order preserved" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_change_tag_dynamic )
{
    using lagrange_t = Lagrange<Dynamic, Scalar>;
    using lagrange_tagged_t = typename lagrange_t::template ChangeTag<2>::type;

    static_assert( lagrange_t::TAG == 0, "Original has TAG 0" );
    static_assert( lagrange_tagged_t::TAG == 2, "Tagged has TAG 2" );
    static_assert( lagrange_tagged_t::is_order_dynamic, "Dynamic preserved" );
}

//=============================================================================
// SECTION 6: Component Basis Tests
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_component_basis_static )
{
    using lagrange_vec_t = Lagrange<2, Vectorial>;
    using component_t = typename lagrange_vec_t::component_basis_type;

    // Component basis should be scalar with same order
    static_assert( component_t::nOrder == 2, "Same order" );

    using fe_comp_t = typename component_t::template apply<2, 2>::type;
    static_assert( fe_comp_t::is_scalar, "Component is scalar" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_component_basis_dynamic )
{
    using lagrange_vec_t = Lagrange<Dynamic, Vectorial>;
    using component_t = typename lagrange_vec_t::component_basis_type;

    // Component basis should be scalar and dynamic
    static_assert( component_t::is_order_dynamic, "Component is dynamic" );

    using fe_comp_t = typename component_t::template apply<2, 2>::type;
    static_assert( fe_comp_t::is_scalar, "Component is scalar" );
}

//=============================================================================
// SECTION 7: Static vs Runtime Numerical Equivalence (P1/P2/P3, 1D/2D/3D)
//=============================================================================

BOOST_AUTO_TEST_CASE( test_lagrange_static_vs_runtime_eval_1d_p1_p2_p3 )
{
    checkStaticRuntimeLagrangeAgreement<1, 1>( 1e-12 );
    checkStaticRuntimeLagrangeAgreement<1, 2>( 1e-11 );
    checkStaticRuntimeLagrangeAgreement<1, 3>( 1e-10 );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_vs_runtime_eval_2d_p1_p2_p3 )
{
    checkStaticRuntimeLagrangeAgreement<2, 1>( 1e-12 );
    checkStaticRuntimeLagrangeAgreement<2, 2>( 1e-11 );
    checkStaticRuntimeLagrangeAgreement<2, 3>( 1e-10 );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_vs_runtime_eval_3d_p1_p2_p3 )
{
    checkStaticRuntimeLagrangeAgreement<3, 1>( 1e-12 );
    checkStaticRuntimeLagrangeAgreement<3, 2>( 1e-11 );
    checkStaticRuntimeLagrangeAgreement<3, 3>( 1e-10 );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_vs_runtime_default_pointset_simplex_high_order )
{
    checkStaticRuntimeLagrangeAgreementDefaultPointSet<2, 4, Simplex<2>>( 1e-11, "default(simplex/Fekete)" );
    checkStaticRuntimeLagrangeAgreementDefaultPointSet<2, 10, Simplex<2>>( 1e-8, "default(simplex/Fekete)" );
}

BOOST_AUTO_TEST_CASE( test_lagrange_static_vs_runtime_default_pointset_hypercube_gausslobatto )
{
    checkStaticRuntimeLagrangeAgreementDefaultPointSet<2, 4, Hypercube<2>>( 1e-12, "default(hypercube/GaussLobatto)" );
    checkStaticRuntimeLagrangeAgreementDefaultPointSet<2, 10, Hypercube<2>>( 1e-10, "default(hypercube/GaussLobatto)" );
}

BOOST_AUTO_TEST_SUITE_END()
