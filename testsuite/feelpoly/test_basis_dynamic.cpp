/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE test_basis_dynamic
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/mp11/algorithm.hpp>

#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/legendre.hpp>

#include <algorithm>
#include <array>
#include <cmath>

using namespace Feel;
namespace bdata = boost::unit_test::data;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Reference dimensions covered by the basis datasets. */
constexpr std::array<unsigned, 3> dimensions{ 1, 2, 3 };
/** Template degrees from which lower/equal/higher-order values are requested. */
constexpr std::array<unsigned, 2> seedOrders{ 2, 3 };
/** Requested degrees used for numerical comparisons against their own bases. */
constexpr std::array<unsigned, 5> requestedOrders{ 0, 1, 2, 3, 4 };
/** Effective degrees used to check runtime table dimensions. */
constexpr std::array<unsigned, 6> runtimeOrders{ 0, 1, 2, 3, 5, 10 };
/** Select simplex Dubiner (`false`) or tensor-product Legendre (`true`). */
constexpr std::array<bool, 2> tensorProducts{ false, true };

/** @brief Assert matching matrix extents, finite entries and bounded maximum error.
 * @tparam A Actual matrix type.
 * @tparam B Reference matrix type.
 * @param actual Table computed through the requested-order overload.
 * @param expected Table computed by the independently instantiated degree.
 */
template<typename A, typename B>
void checkMatrices( A const& actual, B const& expected )
{
    BOOST_REQUIRE_EQUAL( actual.size1(), expected.size1() );
    BOOST_REQUIRE_EQUAL( actual.size2(), expected.size2() );
    double error = 0;
    for ( std::size_t i = 0; i < actual.size1(); ++i )
        for ( std::size_t j = 0; j < actual.size2(); ++j )
        {
            BOOST_REQUIRE( std::isfinite( actual( i, j ) ) );
            error = std::max( error, std::abs( actual( i, j ) - expected( i, j ) ) );
        }
    BOOST_CHECK_SMALL( error, 1e-11 );
}

/** @brief Build three asymmetric interior points for the basis's reference cell.
 * @tparam Basis Basis providing the reference dimension and cell family.
 * @return Matrix with spatial rows and sample columns.
 */
template<typename Basis>
auto samplePoints()
{
    typename Basis::points_type points( Basis::nDim, 3 );
    for ( unsigned d = 0; d < Basis::nDim; ++d )
        for ( unsigned q = 0; q < 3; ++q )
            points( d, q ) = Basis::convex_type::is_simplex
                ? -0.9 + 0.1 * q + 0.025 * d
                : 0.2 + 0.05 * q + 0.1 * d;
    return points;
}

/** @brief Compare requested-order values and every reference derivative direction.
 * @tparam Basis Seed basis whose template degree may differ from the request.
 * @tparam Order Requested degree, instantiated separately for the reference path.
 */
template<typename Basis, unsigned Order>
void checkRequestedOrder()
{
    using requested_type = typename Basis::traits_type::template ChangeOrder<Order>::type;
    auto const points = samplePoints<Basis>();
    checkMatrices( Basis::evaluate( points, Order ), requested_type::evaluate( points ) );
    auto const actual = Basis::derivate( points, Order );
    auto const expected = requested_type::derivate( points );
    BOOST_REQUIRE_EQUAL( actual.size(), static_cast<unsigned>( Basis::nDim ) );
    BOOST_REQUIRE_EQUAL( expected.size(), static_cast<unsigned>( Basis::nDim ) );
    for ( unsigned direction = 0; direction < Basis::nDim; ++direction )
        BOOST_TEST_CONTEXT( "derivative direction=" << direction )
            checkMatrices( actual[direction], expected[direction] );
}

/** @brief Check value/gradient cardinality and point axes at a runtime degree.
 * @tparam Basis Seed basis exposing runtime evaluation overloads.
 * @param order Effective degree whose cardinality determines the basis rows.
 */
template<typename Basis>
void checkShapes( unsigned order )
{
    auto const points = samplePoints<Basis>();
    auto const values = Basis::evaluate( points, order );
    auto const derivatives = Basis::derivate( points, order );
    BOOST_REQUIRE_EQUAL( values.size1(), Basis::convex_type::polyDims( order ) );
    BOOST_REQUIRE_EQUAL( values.size2(), points.size2() );
    BOOST_REQUIRE_EQUAL( derivatives.size(), static_cast<unsigned>( Basis::nDim ) );
    for ( unsigned direction = 0; direction < Basis::nDim; ++direction )
    {
        BOOST_CHECK_EQUAL( derivatives[direction].size1(), values.size1() );
        BOOST_CHECK_EQUAL( derivatives[direction].size2(), values.size2() );
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( basis_dynamic_suite )

/** @test Requested degrees preserve the mode ordering of their own basis. */
BOOST_DATA_TEST_CASE( requested_order_matches_its_own_basis,
    bdata::make( dimensions ) * bdata::make( seedOrders ) *
    bdata::make( requestedOrders ) * bdata::make( tensorProducts ),
    dimension, seedOrder, requestedOrder, tensorProduct )
{
    boost::mp11::mp_with_index<3>( dimension - 1, [&]( auto dimensionTag )
    {
        constexpr unsigned Dim = decltype( dimensionTag )::value + 1;
        boost::mp11::mp_with_index<2>( seedOrder - 2, [&]( auto seedTag )
        {
            constexpr unsigned Seed = decltype( seedTag )::value + 2;
            boost::mp11::mp_with_index<5>( requestedOrder, [&]( auto orderTag )
            {
                constexpr unsigned Order = decltype( orderTag )::value;
                if ( tensorProduct )
                    checkRequestedOrder<Legendre<Dim, Dim, Seed>, Order>();
                else
                    checkRequestedOrder<Dubiner<Dim, Dim, Seed>, Order>();
            } );
        } );
    } );
}

/** @test Runtime value and gradient tables retain their mathematical axes. */
BOOST_DATA_TEST_CASE( runtime_basis_shapes,
    bdata::make( dimensions ) * bdata::make( runtimeOrders ) * bdata::make( tensorProducts ),
    dimension, order, tensorProduct )
{
    boost::mp11::mp_with_index<3>( dimension - 1, [&]( auto dimensionTag )
    {
        constexpr unsigned Dim = decltype( dimensionTag )::value + 1;
        if ( tensorProduct )
            checkShapes<Legendre<Dim, Dim, 3>>( order );
        else
            checkShapes<Dubiner<Dim, Dim, 3>>( order );
    } );
}

BOOST_AUTO_TEST_SUITE_END()
