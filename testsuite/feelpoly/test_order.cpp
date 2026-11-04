/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file test_order.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2026-01-04
   \brief Tests for order system and DOF calculations in feelpoly
 */

#define BOOST_TEST_MODULE test_order
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelmesh/simplex.hpp>
#include <feel/feelpoly/order.hpp>
#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/doftopology.hpp>

#include <boost/mp11/list.hpp>
#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Legacy static provider with no static/runtime order flags. */
struct LegacyStaticOrderProvider
{
    /** Compile-time degree exposed by the legacy contract. */
    static constexpr int nOrder = 2;
};

/** @brief Check substitution-safe detection and validation for a provider type.
 * @tparam Provider Legacy or flagged provider, optionally cv/ref-qualified.
 */
template<typename Provider>
void checkOrderProvider()
{
    using provider_type = std::remove_cvref_t<Provider>;
    constexpr bool dynamic = std::is_same_v<provider_type, OrderBase<Dynamic>>;
    static_assert( hasOrderDynamicFlag<Provider> ==
                   !std::is_same_v<provider_type, LegacyStaticOrderProvider> );
    static_assert( orderIsDynamic<Provider> == dynamic );
    static_assert( orderIsStatic<Provider> == !dynamic );
    BOOST_CHECK_EQUAL( validatedOrder<Provider>( RuntimeOrder{ 2 } ).value, 2 );
    if constexpr ( dynamic )
    {
        static_assert( defaultOrder<Provider>().value == 0 );
        BOOST_CHECK_EQUAL( validatedOrder<Provider>( RuntimeOrder{ 3 } ).value, 3 );
    }
    else
    {
        static_assert( staticOrderValue<Provider>() == 2 );
        static_assert( defaultOrder<Provider>().value == 2 );
        BOOST_CHECK_THROW( (void)validatedOrder<Provider>( RuntimeOrder{ 3 } ), std::invalid_argument );
    }
}

template<typename StaticShape, typename DynamicShape, typename Cardinality>
void checkShapeOrderAndCardinalityOwnership()
{
    static_assert( StaticShape::order_spec >= 0 );
    static_assert( DynamicShape::order_spec == Dynamic );
    static_assert( std::is_base_of_v<OrderBase<StaticShape::order_spec>, StaticShape> );
    static_assert( std::is_base_of_v<OrderBase<Dynamic>, DynamicShape> );
    static_assert(
        StaticShape::numPoints_static ==
        Cardinality::total( StaticShape::nDim, StaticShape::order_spec ) );
    static_assert(
        StaticShape::nbPtsPerVertex_static ==
        Cardinality::perVertex( StaticShape::nDim, StaticShape::order_spec ) );
    static_assert(
        StaticShape::nbPtsPerEdge_static ==
        Cardinality::perEdge( StaticShape::nDim, StaticShape::order_spec ) );
    static_assert(
        StaticShape::nbPtsPerFace_static ==
        Cardinality::perFace( StaticShape::nDim, StaticShape::order_spec ) );
    static_assert(
        StaticShape::nbPtsPerVolume_static ==
        Cardinality::perVolume( StaticShape::nDim, StaticShape::order_spec ) );

    StaticShape staticShape;
    DynamicShape dynamicShape( RuntimeOrder{ StaticShape::order_spec } );
    BOOST_CHECK_EQUAL( staticShape.order(), dynamicShape.order() );
    BOOST_CHECK_EQUAL( staticShape.nPointsTotal(), dynamicShape.nPointsTotal() );
    BOOST_CHECK_EQUAL( staticShape.nPointsOnVertex(), dynamicShape.nPointsOnVertex() );
    BOOST_CHECK_EQUAL( staticShape.nPointsOnEdge(), dynamicShape.nPointsOnEdge() );
    BOOST_CHECK_EQUAL( staticShape.nPointsOnFace(), dynamicShape.nPointsOnFace() );
    BOOST_CHECK_EQUAL( staticShape.nPointsOnVolume(), dynamicShape.nPointsOnVolume() );

    dynamicShape.clearDetachedOrder();
    BOOST_CHECK_EQUAL( dynamicShape.order(), 0 );
}

template<typename StaticTopology, typename DynamicTopology>
void checkTopologyOrderAndCardinalityOwnership()
{
    static_assert( StaticTopology::order_spec >= 0 );
    static_assert( DynamicTopology::order_spec == Dynamic );
    static_assert( std::is_empty_v<StaticTopology> );
    static_assert( sizeof( DynamicTopology ) == sizeof( RuntimeOrder ) );

    constexpr StaticTopology staticTopology;
    constexpr DynamicTopology dynamicTopology(
        RuntimeOrder{ StaticTopology::order_spec } );
    static_assert( staticTopology.localDof() == dynamicTopology.localDof() );
    static_assert( staticTopology.dofPerVertex() == dynamicTopology.dofPerVertex() );
    static_assert( staticTopology.dofPerEdge() == dynamicTopology.dofPerEdge() );
    static_assert( staticTopology.dofPerFace() == dynamicTopology.dofPerFace() );
    static_assert( staticTopology.dofPerVolume() == dynamicTopology.dofPerVolume() );

    BOOST_CHECK_EQUAL( staticTopology.localDof(), dynamicTopology.localDof() );
}
} // namespace

BOOST_AUTO_TEST_SUITE( order_suite )

/** Static legacy, flagged static and flagged dynamic metadata providers. */
using order_providers = boost::mp11::mp_list<LegacyStaticOrderProvider, OrderBase<2>, OrderBase<Dynamic>>;

/** @test Order queries support legacy providers and retain cv/ref-independent flags. */
BOOST_AUTO_TEST_CASE_TEMPLATE( order_provider_detection_is_substitution_safe, Provider, order_providers )
{
    static_assert( !orderIsDynamic<int> );
    checkOrderProvider<Provider>();
    checkOrderProvider<Provider const&>();
    checkOrderProvider<Provider volatile&&>();
}

//
// Test OrderTag basic properties
//
BOOST_AUTO_TEST_CASE( test_order_t_static )
{
    static_assert( std::is_same_v<OrderTag<2>, order_t<2>> );
    static_assert( OrderTag<0>::order_spec == 0 );
    static_assert( OrderTag<5>::order_spec == 5 );
    // Static order values
    static_assert( OrderTag<0>::value == 0, "OrderTag<0>::value should be 0" );
    static_assert( OrderTag<1>::value == 1, "OrderTag<1>::value should be 1" );
    static_assert( OrderTag<5>::value == 5, "OrderTag<5>::value should be 5" );
    static_assert( OrderTag<10>::value == 10, "OrderTag<10>::value should be 10" );

    // Static/dynamic flags for static orders
    static_assert( OrderTag<0>::is_static, "OrderTag<0> should be static" );
    static_assert( OrderTag<3>::is_static, "OrderTag<3> should be static" );
    static_assert( !OrderTag<0>::is_dynamic, "OrderTag<0> should not be dynamic" );
    static_assert( !OrderTag<3>::is_dynamic, "OrderTag<3> should not be dynamic" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_order_t_dynamic )
{
    static_assert( OrderTag<Dynamic>::order_spec == Dynamic );
    // Dynamic order
    static_assert( OrderTag<Dynamic>::value == Dynamic, "OrderTag<Dynamic>::value should be Dynamic" );
    static_assert( OrderTag<Dynamic>::value == -1, "OrderTag<Dynamic>::value should be -1" );

    // Static/dynamic flags for dynamic order
    static_assert( OrderTag<Dynamic>::is_dynamic, "OrderTag<Dynamic> should be dynamic" );
    static_assert( !OrderTag<Dynamic>::is_static, "OrderTag<Dynamic> should not be static" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_runtime_order_checked )
{
    const RuntimeOrder valid = RuntimeOrder::checked( 3 );
    BOOST_CHECK_EQUAL( static_cast<uint16_type>( valid ), 3 );
    BOOST_CHECK_THROW( []() { (void)RuntimeOrder::checked( -1 ); }(), std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( test_order_base )
{
    static_assert( OrderBase<2>::order_spec == 2 );
    static_assert( OrderBase<Dynamic>::order_spec == Dynamic );
    static_assert( std::is_empty_v<OrderBase<2>> );
    static_assert( sizeof( OrderBase<Dynamic> ) == sizeof( RuntimeOrder ) );
    static_assert( OrderBase<2>::is_order_static, "OrderBase<2> should be static" );
    static_assert( !OrderBase<2>::is_order_dynamic, "OrderBase<2> should not be dynamic" );
    static_assert( OrderBase<Dynamic>::is_order_dynamic, "OrderBase<Dynamic> should be dynamic" );
    static_assert( !OrderBase<Dynamic>::is_order_static, "OrderBase<Dynamic> should not be static" );

    constexpr OrderBase<2> staticOrder;
    static_assert( staticOrder.order() == 2, "static order mismatch" );
    static_assert( staticOrder.order() == 2, "static order mismatch" );
    constexpr OrderBase<2> matchingStaticOrder( RuntimeOrder{ 2 } );
    static_assert( matchingStaticOrder.order() == 2 );
    BOOST_CHECK_THROW( ( OrderBase<2>{ RuntimeOrder{ 3 } } ),
                       std::invalid_argument );

    const OrderBase<Dynamic> dynamicOrder( RuntimeOrder( 4 ) );
    BOOST_CHECK_EQUAL( dynamicOrder.order(), 4 );
    BOOST_CHECK_EQUAL( dynamicOrder.order(), 4 );
}

BOOST_AUTO_TEST_CASE( test_authoritative_order_and_cardinality_ownership )
{
    checkShapeOrderAndCardinalityOwnership<
        Simplex<3, 3>, Simplex<3, Dynamic>, detail::SimplexCardinality>();
    checkShapeOrderAndCardinalityOwnership<
        Hypercube<3, 3>, Hypercube<3, Dynamic>, detail::HypercubeCardinality>();

    checkTopologyOrderAndCardinalityOwnership<
        DofTopology<SimplexTag, 3, 3>,
        DofTopology<SimplexTag, 3, Dynamic>>();
    checkTopologyOrderAndCardinalityOwnership<
        DofTopology<HypercubeTag, 3, 3>,
        DofTopology<HypercubeTag, 3, Dynamic>>();
}

BOOST_AUTO_TEST_CASE( test_doftopology_simplex )
{
    using simplex_static_t = DofTopology<SimplexTag, 2, 2>;
    constexpr simplex_static_t simplexStatic;
    static_assert( simplex_static_t::is_order_static );
    static_assert( simplexStatic.order() == 2 );
    static_assert( simplexStatic.localDof() == 6 );
    static_assert( simplexStatic.dofPerVertex() == 1 );
    static_assert( simplexStatic.dofPerEdge() == 1 );
    static_assert( simplexStatic.dofPerFace() == 0 );

    const DofTopology<SimplexTag, 3, Dynamic> simplexDynamic( RuntimeOrder( 3 ) );
    BOOST_CHECK_EQUAL( simplexDynamic.order(), 3 );
    BOOST_CHECK_EQUAL( simplexDynamic.localDof(), 20 );
    BOOST_CHECK_EQUAL( simplexDynamic.dofPerVertex(), 1 );
    BOOST_CHECK_EQUAL( simplexDynamic.dofPerEdge(), 2 );
    BOOST_CHECK_EQUAL( simplexDynamic.dofPerFace(), 1 );
    BOOST_CHECK_EQUAL( simplexDynamic.dofPerVolume(), 0 );
    BOOST_CHECK_EQUAL( simplexDynamic.dofPerEntity( 1, 0 ), simplexDynamic.dofPerEdge() );
}

BOOST_AUTO_TEST_CASE( test_doftopology_hypercube )
{
    using hypercube_static_t = DofTopology<HypercubeTag, 2, 2>;
    constexpr hypercube_static_t hypercubeStatic;
    static_assert( hypercube_static_t::is_order_static );
    static_assert( hypercubeStatic.order() == 2 );
    static_assert( hypercubeStatic.localDof() == 9 );
    static_assert( hypercubeStatic.dofPerVertex() == 1 );
    static_assert( hypercubeStatic.dofPerEdge() == 1 );
    static_assert( hypercubeStatic.dofPerFace() == 1 );

    const DofTopology<HypercubeTag, 3, Dynamic> hypercubeDynamic( RuntimeOrder( 2 ) );
    BOOST_CHECK_EQUAL( hypercubeDynamic.order(), 2 );
    BOOST_CHECK_EQUAL( hypercubeDynamic.localDof(), 27 );
    BOOST_CHECK_EQUAL( hypercubeDynamic.dofPerVertex(), 1 );
    BOOST_CHECK_EQUAL( hypercubeDynamic.dofPerEdge(), 1 );
    BOOST_CHECK_EQUAL( hypercubeDynamic.dofPerFace(), 1 );
    BOOST_CHECK_EQUAL( hypercubeDynamic.dofPerVolume(), 1 );
    BOOST_CHECK_EQUAL( hypercubeDynamic.dofPerEntity( 3, 0 ), hypercubeDynamic.dofPerVolume() );
}

//
// Test binomial coefficient calculation
//
BOOST_AUTO_TEST_CASE( test_binomial )
{
    // Known binomial coefficients C(n, k) = n! / (k! * (n-k)!)
    // C(0, 0) = 1
    static_assert( detail::binomial( 0, 0 ) == 1, "C(0,0) = 1" );

    // C(n, 0) = 1 for all n
    static_assert( detail::binomial( 1, 0 ) == 1, "C(1,0) = 1" );
    static_assert( detail::binomial( 5, 0 ) == 1, "C(5,0) = 1" );
    static_assert( detail::binomial( 10, 0 ) == 1, "C(10,0) = 1" );

    // C(n, n) = 1 for all n
    static_assert( detail::binomial( 1, 1 ) == 1, "C(1,1) = 1" );
    static_assert( detail::binomial( 5, 5 ) == 1, "C(5,5) = 1" );
    static_assert( detail::binomial( 10, 10 ) == 1, "C(10,10) = 1" );

    // C(n, 1) = n
    static_assert( detail::binomial( 5, 1 ) == 5, "C(5,1) = 5" );
    static_assert( detail::binomial( 10, 1 ) == 10, "C(10,1) = 10" );

    // Pascal's triangle values
    static_assert( detail::binomial( 4, 2 ) == 6, "C(4,2) = 6" );
    static_assert( detail::binomial( 5, 2 ) == 10, "C(5,2) = 10" );
    static_assert( detail::binomial( 6, 3 ) == 20, "C(6,3) = 20" );
    static_assert( detail::binomial( 10, 5 ) == 252, "C(10,5) = 252" );

    // Invalid inputs
    static_assert( detail::binomial( 3, 5 ) == 0, "C(3,5) = 0 (k > n)" );
    static_assert( detail::binomial( 5, -1 ) == 0, "C(5,-1) = 0 (k < 0)" );

    BOOST_CHECK( true );
}

//
// Test pow_int calculation
//
BOOST_AUTO_TEST_CASE( test_pow_int )
{
    // base^0 = 1
    static_assert( detail::powInt( 2, 0 ) == 1, "2^0 = 1" );
    static_assert( detail::powInt( 5, 0 ) == 1, "5^0 = 1" );

    // base^1 = base
    static_assert( detail::powInt( 2, 1 ) == 2, "2^1 = 2" );
    static_assert( detail::powInt( 5, 1 ) == 5, "5^1 = 5" );

    // Powers of 2
    static_assert( detail::powInt( 2, 2 ) == 4, "2^2 = 4" );
    static_assert( detail::powInt( 2, 3 ) == 8, "2^3 = 8" );
    static_assert( detail::powInt( 2, 4 ) == 16, "2^4 = 16" );
    static_assert( detail::powInt( 2, 10 ) == 1024, "2^10 = 1024" );

    // Other bases
    static_assert( detail::powInt( 3, 3 ) == 27, "3^3 = 27" );
    static_assert( detail::powInt( 4, 3 ) == 64, "4^3 = 64" );
    static_assert( detail::powInt( 5, 4 ) == 625, "5^4 = 625" );

    BOOST_CHECK( true );
}

//
// Test simplex DOF calculations
// DOF for simplex of dimension d with order p: C(d+p, d) = (d+p)! / (d! * p!)
//
BOOST_AUTO_TEST_CASE( test_simplex_dof_static )
{
    // 1D simplex (segment)
    // P0: 1 DOF, P1: 2 DOFs, P2: 3 DOFs, Pk: k+1 DOFs
    static_assert( detail::simplexDofStatic<1, 0>() == 1, "1D P0: 1 DOF" );
    static_assert( detail::simplexDofStatic<1, 1>() == 2, "1D P1: 2 DOFs" );
    static_assert( detail::simplexDofStatic<1, 2>() == 3, "1D P2: 3 DOFs" );
    static_assert( detail::simplexDofStatic<1, 3>() == 4, "1D P3: 4 DOFs" );
    static_assert( detail::simplexDofStatic<1, 5>() == 6, "1D P5: 6 DOFs" );

    // 2D simplex (triangle)
    // P0: 1, P1: 3, P2: 6, P3: 10, Pk: (k+1)(k+2)/2
    static_assert( detail::simplexDofStatic<2, 0>() == 1, "2D P0: 1 DOF" );
    static_assert( detail::simplexDofStatic<2, 1>() == 3, "2D P1: 3 DOFs" );
    static_assert( detail::simplexDofStatic<2, 2>() == 6, "2D P2: 6 DOFs" );
    static_assert( detail::simplexDofStatic<2, 3>() == 10, "2D P3: 10 DOFs" );
    static_assert( detail::simplexDofStatic<2, 4>() == 15, "2D P4: 15 DOFs" );
    static_assert( detail::simplexDofStatic<2, 5>() == 21, "2D P5: 21 DOFs" );

    // 3D simplex (tetrahedron)
    // P0: 1, P1: 4, P2: 10, P3: 20, Pk: (k+1)(k+2)(k+3)/6
    static_assert( detail::simplexDofStatic<3, 0>() == 1, "3D P0: 1 DOF" );
    static_assert( detail::simplexDofStatic<3, 1>() == 4, "3D P1: 4 DOFs" );
    static_assert( detail::simplexDofStatic<3, 2>() == 10, "3D P2: 10 DOFs" );
    static_assert( detail::simplexDofStatic<3, 3>() == 20, "3D P3: 20 DOFs" );
    static_assert( detail::simplexDofStatic<3, 4>() == 35, "3D P4: 35 DOFs" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_simplex_dof_dynamic )
{
    // Compare dynamic calculation with static values
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<1>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<1>( 1 ), 2 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<1>( 2 ), 3 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<1>( 5 ), 6 );

    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<2>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<2>( 1 ), 3 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<2>( 2 ), 6 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<2>( 3 ), 10 );

    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<3>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<3>( 1 ), 4 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<3>( 2 ), 10 );
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<3>( 3 ), 20 );

    // Invalid order
    BOOST_CHECK_EQUAL( detail::simplexDofDynamic<2>( -1 ), 0 );
}

//
// Test hypercube DOF calculations
// DOF for hypercube of dimension d with order p: (p+1)^d
//
BOOST_AUTO_TEST_CASE( test_hypercube_dof_static )
{
    // 1D hypercube (segment) - same as simplex
    static_assert( detail::hypercubeDofStatic<1, 0>() == 1, "1D Q0: 1 DOF" );
    static_assert( detail::hypercubeDofStatic<1, 1>() == 2, "1D Q1: 2 DOFs" );
    static_assert( detail::hypercubeDofStatic<1, 2>() == 3, "1D Q2: 3 DOFs" );
    static_assert( detail::hypercubeDofStatic<1, 3>() == 4, "1D Q3: 4 DOFs" );

    // 2D hypercube (quadrilateral)
    // Q0: 1, Q1: 4, Q2: 9, Q3: 16, Qk: (k+1)^2
    static_assert( detail::hypercubeDofStatic<2, 0>() == 1, "2D Q0: 1 DOF" );
    static_assert( detail::hypercubeDofStatic<2, 1>() == 4, "2D Q1: 4 DOFs" );
    static_assert( detail::hypercubeDofStatic<2, 2>() == 9, "2D Q2: 9 DOFs" );
    static_assert( detail::hypercubeDofStatic<2, 3>() == 16, "2D Q3: 16 DOFs" );
    static_assert( detail::hypercubeDofStatic<2, 4>() == 25, "2D Q4: 25 DOFs" );

    // 3D hypercube (hexahedron)
    // Q0: 1, Q1: 8, Q2: 27, Q3: 64, Qk: (k+1)^3
    static_assert( detail::hypercubeDofStatic<3, 0>() == 1, "3D Q0: 1 DOF" );
    static_assert( detail::hypercubeDofStatic<3, 1>() == 8, "3D Q1: 8 DOFs" );
    static_assert( detail::hypercubeDofStatic<3, 2>() == 27, "3D Q2: 27 DOFs" );
    static_assert( detail::hypercubeDofStatic<3, 3>() == 64, "3D Q3: 64 DOFs" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_hypercube_dof_dynamic )
{
    // Compare dynamic calculation with static values
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<1>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<1>( 1 ), 2 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<1>( 2 ), 3 );

    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<2>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<2>( 1 ), 4 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<2>( 2 ), 9 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<2>( 3 ), 16 );

    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<3>( 0 ), 1 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<3>( 1 ), 8 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<3>( 2 ), 27 );
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<3>( 3 ), 64 );

    // Invalid order
    BOOST_CHECK_EQUAL( detail::hypercubeDofDynamic<2>( -1 ), 0 );
}

//
// Test DofCalculator with static order
//
BOOST_AUTO_TEST_CASE( test_dof_calculator_static )
{
    // 2D simplex, order 3
    using calc_2d_3 = DofCalculator<2, OrderTag<3>>;
    static_assert( calc_2d_3::simplex() == 10, "2D P3 simplex: 10 DOFs" );
    static_assert( calc_2d_3::hypercube() == 16, "2D Q3 hypercube: 16 DOFs" );

    // 3D simplex, order 2
    using calc_3d_2 = DofCalculator<3, OrderTag<2>>;
    static_assert( calc_3d_2::simplex() == 10, "3D P2 simplex: 10 DOFs" );
    static_assert( calc_3d_2::hypercube() == 27, "3D Q2 hypercube: 27 DOFs" );

    // 1D, order 5
    using calc_1d_5 = DofCalculator<1, OrderTag<5>>;
    static_assert( calc_1d_5::simplex() == 6, "1D P5: 6 DOFs" );
    static_assert( calc_1d_5::hypercube() == 6, "1D Q5: 6 DOFs" );

    BOOST_CHECK( true );
}

//
// Test DofCalculator with dynamic order
//
BOOST_AUTO_TEST_CASE( test_dof_calculator_dynamic )
{
    using calc_2d_dyn = DofCalculator<2, OrderTag<Dynamic>>;
    using calc_3d_dyn = DofCalculator<3, OrderTag<Dynamic>>;

    // 2D tests
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 1 ), 3 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 2 ), 6 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::simplex( 3 ), 10 );

    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 1 ), 4 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 2 ), 9 );
    BOOST_CHECK_EQUAL( calc_2d_dyn::hypercube( 3 ), 16 );

    // 3D tests
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 1 ), 4 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 2 ), 10 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::simplex( 3 ), 20 );

    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 0 ), 1 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 1 ), 8 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 2 ), 27 );
    BOOST_CHECK_EQUAL( calc_3d_dyn::hypercube( 3 ), 64 );
}

//
// Test order_t satisfies concepts
//
BOOST_AUTO_TEST_CASE( test_order_concepts_satisfaction )
{
    // Static orders should satisfy StaticOrder concept
    static_assert( StaticOrder<OrderTag<0>>, "OrderTag<0> should satisfy StaticOrder" );
    static_assert( StaticOrder<OrderTag<1>>, "OrderTag<1> should satisfy StaticOrder" );
    static_assert( StaticOrder<OrderTag<5>>, "OrderTag<5> should satisfy StaticOrder" );

    // Dynamic order should satisfy DynamicOrder concept
    static_assert( DynamicOrder<OrderTag<Dynamic>>, "OrderTag<Dynamic> should satisfy DynamicOrder" );

    // Cross-checks
    static_assert( !DynamicOrder<OrderTag<0>>, "OrderTag<0> should NOT satisfy DynamicOrder" );
    static_assert( !StaticOrder<OrderTag<Dynamic>>, "OrderTag<Dynamic> should NOT satisfy StaticOrder" );

    // All should satisfy HasOrder
    static_assert( HasOrder<OrderTag<0>>, "OrderTag<0> should satisfy HasOrder" );
    static_assert( HasOrder<OrderTag<5>>, "OrderTag<5> should satisfy HasOrder" );
    static_assert( HasOrder<OrderTag<Dynamic>>, "OrderTag<Dynamic> should satisfy HasOrder" );

    BOOST_CHECK( true );
}

//
// Test comparison of simplex vs hypercube DOFs
//
BOOST_AUTO_TEST_CASE( test_simplex_vs_hypercube_dofs )
{
    // For same dimension and order, hypercube always has >= DOFs than simplex
    // (except in 1D where they are equal)

    // 1D: equal
    static_assert( detail::simplexDofStatic<1, 3>() == detail::hypercubeDofStatic<1, 3>(),
                   "1D: simplex and hypercube have same DOFs" );

    // 2D: hypercube has more
    static_assert( detail::hypercubeDofStatic<2, 2>() > detail::simplexDofStatic<2, 2>(),
                   "2D: hypercube has more DOFs (9 > 6)" );
    static_assert( detail::hypercubeDofStatic<2, 3>() > detail::simplexDofStatic<2, 3>(),
                   "2D: hypercube has more DOFs (16 > 10)" );

    // 3D: hypercube has more
    static_assert( detail::hypercubeDofStatic<3, 2>() > detail::simplexDofStatic<3, 2>(),
                   "3D: hypercube has more DOFs (27 > 10)" );
    static_assert( detail::hypercubeDofStatic<3, 3>() > detail::simplexDofStatic<3, 3>(),
                   "3D: hypercube has more DOFs (64 > 20)" );

    BOOST_CHECK( true );
}

//
// Test high-order DOF calculations
//
BOOST_AUTO_TEST_CASE( test_high_order_dofs )
{
    // Test that high-order calculations don't overflow for reasonable orders
    // 2D simplex P10: C(12, 2) = 66
    static_assert( detail::simplexDofStatic<2, 10>() == 66, "2D P10: 66 DOFs" );

    // 3D simplex P5: C(8, 3) = 56
    static_assert( detail::simplexDofStatic<3, 5>() == 56, "3D P5: 56 DOFs" );

    // 2D hypercube Q10: 11^2 = 121
    static_assert( detail::hypercubeDofStatic<2, 10>() == 121, "2D Q10: 121 DOFs" );

    // 3D hypercube Q5: 6^3 = 216
    static_assert( detail::hypercubeDofStatic<3, 5>() == 216, "3D Q5: 216 DOFs" );

    BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_local_index_cardinality_boundaries )
{
    using simplex_topology = DofTopology<SimplexTag, 3, Dynamic>;
    using hypercube_topology = DofTopology<HypercubeTag, 3, Dynamic>;
    using simplex_shape = Simplex<3, Dynamic>;
    using hypercube_shape = Hypercube<3, Dynamic>;

    const simplex_topology simplexAtLimit( RuntimeOrder{ 71 } );
    BOOST_CHECK_EQUAL( simplexAtLimit.localDof(), 64824 );
    BOOST_CHECK_THROW( simplex_topology( RuntimeOrder{ 72 } ),
                       std::length_error );
    BOOST_CHECK_THROW( simplex_shape( RuntimeOrder{ 72 } ),
                       std::length_error );

    const hypercube_topology hypercubeAtLimit( RuntimeOrder{ 39 } );
    BOOST_CHECK_EQUAL( hypercubeAtLimit.localDof(), 64000 );
    BOOST_CHECK_THROW( hypercube_topology( RuntimeOrder{ 40 } ),
                       std::length_error );
    BOOST_CHECK_THROW( hypercube_shape( RuntimeOrder{ 40 } ),
                       std::length_error );

    BOOST_CHECK_EQUAL(
        detail::checkedLocalIndexValue(
            3 * detail::SimplexCardinality::total( 3, 48 ) ),
        62475 );
    BOOST_CHECK_THROW(
        (void)detail::checkedLocalIndexValue(
            3 * detail::SimplexCardinality::total( 3, 49 ) ),
        std::length_error );

    BOOST_CHECK_EQUAL( detail::maxLocalIndexCardinality, 65535 );
    BOOST_CHECK_THROW(
        (void)detail::checkedLocalIndexValue( 65536 ), std::length_error );
}

BOOST_AUTO_TEST_SUITE_END()
