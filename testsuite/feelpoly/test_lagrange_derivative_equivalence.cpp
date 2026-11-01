/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file test_lagrange_derivative_equivalence.cpp
   \brief Dedicated tests for compile-time vs runtime Lagrange derivative equivalence
 */

#define BOOST_TEST_MODULE test_lagrange_derivative_equivalence
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelmesh/simplex.hpp>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <string>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
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

template <typename FEStaticType, typename FEDynamicType, typename PointsType>
void checkDerivativeAgreement( FEStaticType const& fe_static,
                               FEDynamicType const& fe_dynamic,
                               PointsType const& pts,
                               double tol,
                               std::string const& label )
{
    auto const d_static = fe_static.derivate( pts );
    auto const d_dynamic = fe_dynamic.derivate( pts );

    BOOST_TEST_CONTEXT( label + " derivative vector size" )
    {
        BOOST_REQUIRE_EQUAL( d_static.size(), d_dynamic.size() );
    }

    for ( std::size_t c = 0; c < d_static.size(); ++c )
    {
        checkMatricesNear( d_static[c],
                           d_dynamic[c],
                           tol,
                           label + " d/dx" + std::to_string( c ) + " via derivate(points)" );

        auto const derivative_static = fe_static.derivate( c );
        auto const derivative_dynamic = fe_dynamic.derivate( c );
        BOOST_TEST_CONTEXT( label + " derivative expansion order" )
        {
            BOOST_CHECK_EQUAL( derivative_static.order(), fe_static.order() );
            BOOST_CHECK_EQUAL( derivative_dynamic.order(), fe_dynamic.order() );
        }
        auto const dpoly_static = derivative_static.evaluate( pts );
        checkMatricesNear( d_static[c],
                           dpoly_static,
                           tol,
                           label + " d/dx" + std::to_string( c ) + " static internal consistency" );

        auto const dpoly_dynamic = derivative_dynamic.evaluate( pts );
        checkMatricesNear( dpoly_static,
                           dpoly_dynamic,
                           tol,
                           label + " d/dx" + std::to_string( c ) + " via derivate(i).evaluate(points)" );
        checkMatricesNear( d_dynamic[c],
                           dpoly_dynamic,
                           tol,
                           label + " d/dx" + std::to_string( c ) + " runtime internal consistency" );
        checkMatricesNear( derivative_dynamic.derivate( c ).evaluate( pts ),
                           fe_dynamic.derivate( c, c, pts ),
                           tol,
                           label + " d2/dx" + std::to_string( c ) + "2 chained runtime derivative" );
    }
}

template <int Dim, int Order>
void checkStaticRuntimeDerivativeAgreement( double tol )
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

    checkDerivativeAgreement( fe_static, fe_dynamic, pts_static, tol, prefix + "derivatives on static nodes" );
    checkDerivativeAgreement( fe_static, fe_dynamic, pts_dynamic, tol, prefix + "derivatives on dynamic nodes" );

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
    checkDerivativeAgreement( fe_static, fe_dynamic, probePts, tol, prefix + "derivatives on probe points" );
}
}

BOOST_AUTO_TEST_SUITE( lagrange_derivative_equivalence_suite )

BOOST_AUTO_TEST_CASE( test_lagrange_derivative_equivalence_1d_p1_p2_p3 )
{
    checkStaticRuntimeDerivativeAgreement<1, 1>( 1e-10 );
    checkStaticRuntimeDerivativeAgreement<1, 2>( 1e-9 );
    checkStaticRuntimeDerivativeAgreement<1, 3>( 1e-8 );
}

BOOST_AUTO_TEST_CASE( test_lagrange_derivative_equivalence_2d_p1_p2_p3 )
{
    checkStaticRuntimeDerivativeAgreement<2, 1>( 1e-10 );
    checkStaticRuntimeDerivativeAgreement<2, 2>( 1e-9 );
    checkStaticRuntimeDerivativeAgreement<2, 3>( 1e-8 );
}

BOOST_AUTO_TEST_CASE( test_lagrange_derivative_equivalence_3d_p1_p2_p3 )
{
    checkStaticRuntimeDerivativeAgreement<3, 1>( 1e-10 );
    checkStaticRuntimeDerivativeAgreement<3, 2>( 1e-9 );
    checkStaticRuntimeDerivativeAgreement<3, 3>( 1e-8 );
}

BOOST_AUTO_TEST_SUITE_END()
