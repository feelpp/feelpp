/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file test_legendre_tabulation.cpp
 * @brief Tensor values, gradients and Hessians against independent Jacobi jets.
 */
#define BOOST_TEST_MODULE legendre_tabulation
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/legendretabulation.hpp>
#include <boost/mp11/list.hpp>
#include "polynomial_reference.hpp"
#include <array>
#include <cmath>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Share independent modal-jet checks across dimensions, degrees and normalization. */
template<int Dim, int P, bool Normalized>
void checkLegendreTabulation()
{
    using Static = TensorLegendreBasis<Dim, P, double, Normalized>;
    using Runtime = TensorLegendreBasis<Dim, Dynamic, double, Normalized>;
    Static fixed;
    Runtime runtime{ RuntimeOrder{ P } };
    Eigen::Matrix<double, Dim, Dim + 3> points;
    points.setConstant( -1. );
    for ( int d = 0; d < Dim; ++d ) points( d, d + 1 ) = 1.;
    points.col( Dim + 1 ).setOnes();
    points.col( Dim + 2 ).setConstant( .25 );
    auto values = fixed.values( points );
    auto dynamicValues = runtime.values( points );
    auto gradients = fixed.gradients( points );
    auto dynamicGradients = runtime.gradients( points );
    auto hessians = fixed.hessians( points );
    auto dynamicHessians = runtime.hessians( points );
    static_assert( decltype( values )::rowsAtCompileTime == Static::basisExtent );
    static_assert( decltype( values )::colsAtCompileTime == Dim + 3 );
    static_assert( decltype( values )::usesInlineStorage );
    static_assert( decltype( dynamicValues )::rowsAtCompileTime == Dynamic );
    static_assert( decltype( dynamicValues )::colsAtCompileTime == Dim + 3 );
    static_assert( !decltype( dynamicValues )::usesInlineStorage );
    BOOST_CHECK_SMALL( ( values.matrix() - dynamicValues.matrix() ).norm(), 1e-11 );
    for ( int d = 0; d < Dim; ++d )
        BOOST_CHECK_SMALL( ( gradients[d].matrix() - dynamicGradients[d].matrix() ).norm(), 1e-11 );
    for ( int ij = 0; ij < Dim * Dim; ++ij )
        BOOST_CHECK_SMALL( ( hessians[ij].matrix() - dynamicHessians[ij].matrix() ).norm(), 1e-10 );
    for ( Eigen::Index row = 0; row < fixed.basisCount(); ++row )
    {
        std::array<unsigned, Dim> degrees;
        auto remaining = row;
        double scale = 1.;
        for ( int d = Dim - 1; d >= 0; --d )
        {
            degrees[d] = remaining % ( P + 1 );
            remaining /= P + 1;
            if constexpr ( Normalized ) scale *= degrees[d] + .5;
        }
        if constexpr ( Normalized ) scale = std::sqrt( scale );
        for ( int q = 0; q < points.cols(); ++q )
        {
            using Jet = decltype( polynomial_test::jacobi( 0, 0., 0., 0. ) );
            std::array<Jet, Dim> factors;
            for ( int d = 0; d < Dim; ++d )
                factors[d] = polynomial_test::jacobi( degrees[d], 0., 0., points( d, q ) );
            auto expected = [&]( int first, int second )
            {
                polynomial_test::precise_type product = scale;
                for ( int d = 0; d < Dim; ++d )
                {
                    int order = ( first == d ) + ( second == d );
                    product *= order == 0 ? factors[d].value :
                               order == 1 ? factors[d].gradient[0] : factors[d].hessian[0];
                }
                return product.template convert_to<double>();
            };
            BOOST_CHECK_SMALL( values.matrix()( row, q ) - expected( -1, -1 ), 1e-10 );
            for ( int i = 0; i < Dim; ++i )
            {
                BOOST_CHECK_SMALL( gradients[i].matrix()( row, q ) - expected( i, -1 ), 1e-9 );
                for ( int j = 0; j < Dim; ++j )
                    BOOST_CHECK_SMALL( hessians[i * Dim + j].matrix()( row, q ) - expected( i, j ), 1e-8 );
            }
        }
    }
    typename Static::template Workspace<0> workspace;
    auto const* address = values.matrix().data();
    fixed.valuesInto( values, points, workspace );
    fixed.valuesInto( values, points, workspace );
    BOOST_CHECK( values.matrix().data() == address );
    if constexpr ( Static::basisExtent >= Dim )
    {
        auto checkAlias = [&]( auto& output, auto& coordinateOwner, auto& work, auto&& evaluate )
        {
            auto expected = output;
            coordinateOwner.matrix().topRows( Dim ) = points;
            evaluate( output, coordinateOwner.matrix().topRows( Dim ), work );
            if constexpr ( requires { output.matrix(); } )
                BOOST_CHECK_SMALL( ( output.matrix() - expected.matrix() ).norm(), 1e-11 );
            else
                for ( std::size_t d = 0; d < output.size(); ++d )
                    BOOST_CHECK_SMALL( ( output[d].matrix() - expected[d].matrix() ).norm(), 1e-10 );
        };
        checkAlias( values, values, workspace, [&]( auto& output, auto const& coordinates, auto& work )
        { fixed.valuesInto( output, coordinates, work ); } );
        typename Static::template Workspace<1> gradientWorkspace;
        checkAlias( gradients, gradients.back(), gradientWorkspace,
            [&]( auto& output, auto const& coordinates, auto& work )
            { fixed.gradientsInto( output, coordinates, work ); } );
        typename Static::template Workspace<2> hessianWorkspace;
        checkAlias( hessians, hessians.back(), hessianWorkspace,
            [&]( auto& output, auto const& coordinates, auto& work )
            { fixed.hessiansInto( output, coordinates, work ); } );
    }
    auto expressionValues = fixed.values( points + points * 0. );
    auto expressionGradients = fixed.gradients( points + points * 0. );
    auto expressionHessians = fixed.hessians( points + points * 0. );
    BOOST_CHECK_SMALL( ( expressionValues.matrix() - values.matrix() ).norm(), 1e-11 );
    for ( int d = 0; d < Dim; ++d )
        BOOST_CHECK_SMALL( ( expressionGradients[d].matrix() - gradients[d].matrix() ).norm(), 1e-10 );
    for ( int ij = 0; ij < Dim * Dim; ++ij )
        BOOST_CHECK_SMALL( ( expressionHessians[ij].matrix() - hessians[ij].matrix() ).norm(), 1e-10 );
    Eigen::Matrix<double, Dim, Dynamic> singlePoint = points.leftCols( 1 );
    auto singleStatic = fixed.values( singlePoint );
    auto singleRuntime = runtime.values( singlePoint );
    static_assert( decltype( singleStatic )::colsAtCompileTime == Dynamic );
    static_assert( decltype( singleRuntime )::colsAtCompileTime == Dynamic );
    static_assert( decltype( singleRuntime )::rowsAtCompileTime == Dynamic );
    BOOST_CHECK_SMALL( ( singleStatic.matrix().col( 0 ) - values.matrix().col( 0 ) ).norm(), 1e-11 );
    BOOST_CHECK_SMALL( ( singleRuntime.matrix().col( 0 ) - values.matrix().col( 0 ) ).norm(), 1e-11 );
    typename Runtime::template Workspace<0> wrong{ RuntimeOrder{ P + 1 } };
    BOOST_CHECK_THROW( runtime.valuesInto( dynamicValues, points, wrong ), std::invalid_argument );
}
} // namespace

/** Fixed reference dimensions for the shared mathematical checks. */
using Dimensions = boost::mp11::mp_list_c<int, 1, 2, 3>;

/** @test Static/runtime tensor jets preserve normalization, ordering and derivatives. */
BOOST_AUTO_TEST_CASE_TEMPLATE( tensorLegendreTabulation, Dimension, Dimensions )
{
    constexpr int dim = Dimension::value;
    checkLegendreTabulation<dim, 0, false>();
    checkLegendreTabulation<dim, 1, false>();
    checkLegendreTabulation<dim, 2, false>();
    checkLegendreTabulation<dim, 3, false>();
    checkLegendreTabulation<dim, 5, false>();
    checkLegendreTabulation<dim, 0, true>();
    checkLegendreTabulation<dim, 1, true>();
    checkLegendreTabulation<dim, 2, true>();
    checkLegendreTabulation<dim, 3, true>();
    checkLegendreTabulation<dim, 5, true>();
}
