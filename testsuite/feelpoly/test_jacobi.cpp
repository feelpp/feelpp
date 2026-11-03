/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2005-2006 EPFL
 * SPDX-FileCopyrightText: 2009 Université de Grenoble 1 (Joseph Fourier)
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define BOOST_TEST_MODULE test_jacobi
#include <feel/feelcore/testsuite.hpp>
#include <boost/test/data/test_case.hpp>

#include <feel/feelpoly/jacobi.hpp>
#include <feel/feelpoly/dubiner.hpp>
#include <feel/feelpoly/orthonormalpolynomialset.hpp>

#include "polynomial_reference.hpp"

#include <array>
#include <cmath>

using namespace Feel;
namespace bdata = boost::unit_test::data;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Scalar and batch Jacobi degrees, including constants and higher orders. */
constexpr std::array<unsigned, 6> degrees{ 0, 1, 2, 3, 5, 10 };
/** First Jacobi parameters, paired elementwise with beta. */
constexpr std::array<double, 3> alpha{ 0.0, 0.25, 1.0 };
/** Second Jacobi parameters, paired elementwise with alpha. */
constexpr std::array<double, 3> beta{ 0.0, 0.75, 2.0 };
/** Reference coordinates including both endpoints and asymmetric interior points. */
constexpr std::array<double, 7> samples{ -1.0, -0.7, -0.2, 0.0, 0.3, 0.9, 1.0 };
/** Empty, single-point and multi-point batch sizes. */
constexpr std::array<unsigned, 3> batchSizes{ 0, 1, samples.size() };
/** Reference dimensions used for simplex value and derivative checks. */
constexpr std::array<unsigned, 3> dimensions{ 1, 2, 3 };
/** Bounded total-degree matrix for the ordinary Dubiner reference comparisons. */
constexpr std::array<unsigned, 5> simplexDegrees{ 0, 1, 2, 3, 5 };
/** Non-normalized and orthonormal modal conventions. */
constexpr std::array<bool, 2> normalization{ false, true };

/** @brief Require a finite value within a reference-scaled absolute tolerance.
 * @param actual Value computed by the production polynomial API.
 * @param expected Independently computed high-precision reference.
 * @param tolerance Absolute tolerance multiplied by `1 + abs(expected)`.
 */
void checkNear( double actual, polynomial_test::precise_type const& expected,
                double tolerance = 2e-11 )
{
    double const reference = expected.convert_to<double>();
    BOOST_REQUIRE( std::isfinite( actual ) );
    BOOST_CHECK_SMALL( actual - reference, tolerance * ( 1.0 + std::abs( reference ) ) );
}

/** @brief Check a scalar Jacobi polynomial and its first derivative.
 * @tparam Polynomial Static- or runtime-order scalar polynomial type.
 * @param polynomial Production polynomial callable.
 * @param degree Effective polynomial degree.
 * @param a First Jacobi parameter.
 * @param b Second Jacobi parameter.
 * @param x Reference coordinate.
 */
template<typename Polynomial>
void checkScalar( Polynomial const& polynomial, unsigned degree,
                  double a, double b, double x )
{
    auto const reference = polynomial_test::jacobi( degree, a, b, x );
    checkNear( polynomial( x ), reference.value );
    checkNear( polynomial.derivate( x ), reference.gradient[0] );
}

/** @brief Compare Dubiner values/gradients at vertices and collapsed-limit points.
 * @tparam Dim Reference dimension, restricted to 1, 2 or 3.
 * @tparam NormalizedBasis Whether both paths use the orthonormal convention.
 * @param degree Effective maximum total degree requested from the shared API.
 */
template<int Dim, bool NormalizedBasis>
void checkDubiner( unsigned degree )
{
    using basis_type = Dubiner<Dim, Dim, 1, Normalized<NormalizedBasis>>;
    // Every vertex, the barycenter, and an approach to the collapsed vertex.
    typename basis_type::points_type points( Dim, Dim + 3 );
    for ( int q = 0; q < Dim + 3; ++q )
        for ( int d = 0; d < Dim; ++d )
            points( d, q ) = -1.0;
    for ( int d = 0; d < Dim; ++d )
    {
        points( d, d + 1 ) = 1.0;
        points( d, Dim + 1 ) = -1.0 + 2.0 / ( Dim + 1 );
        points( d, Dim + 2 ) = -1.0 + 1e-8;
    }
    points( Dim - 1, Dim + 2 ) = 1.0 - Dim * 1e-8;
    auto const values = basis_type::evaluate( points, degree );
    auto const gradients = basis_type::derivate( points, degree );
    BOOST_REQUIRE_EQUAL( gradients.size(), Dim );
    BOOST_REQUIRE_EQUAL( values.size2(), points.size2() );
    for ( int q = 0; q < Dim + 3; ++q )
    {
        std::array<double, Dim> point;
        for ( int d = 0; d < Dim; ++d )
            point[d] = points( d, q );
        auto const expected = polynomial_test::dubiner<Dim>( degree, point, NormalizedBasis );
        BOOST_REQUIRE_EQUAL( values.size1(), expected.size() );
        BOOST_TEST_CONTEXT( "dimension=" << Dim << " degree=" << degree <<
                            " normalized=" << NormalizedBasis << " sample=" << q )
        {
            for ( std::size_t i = 0; i < expected.size(); ++i )
            {
                checkNear( values( i, q ), expected[i].value );
                for ( int d = 0; d < Dim; ++d )
                {
                    BOOST_REQUIRE_EQUAL( gradients[d].size1(), expected.size() );
                    BOOST_REQUIRE_EQUAL( gradients[d].size2(), points.size2() );
                    checkNear( gradients[d]( i, q ), expected[i].gradient[d], 2e-8 );
                }
            }
        }
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE( jacobi )

/** @test Scalar/batch Jacobi values and derivatives match the independent oracle. */
BOOST_DATA_TEST_CASE( batch_and_runtime_jacobi,
    bdata::make( degrees ) * ( bdata::make( alpha ) ^ bdata::make( beta ) ) * bdata::make( batchSizes ),
    degree, a, b, batchSize )
{
    ublas::vector<double> points( batchSize );
    for ( unsigned q = 0; q < batchSize; ++q )
        points[q] = samples[q];
    auto const values = JacobiBatchEvaluation<double>( degree, a, b, points );
    auto const gradients = JacobiBatchDerivation<double>( degree, a, b, points );
    auto const runtimeValues = dyna::JacobiBatchEvaluation( degree, a, b, points );
    auto const runtimeGradients = dyna::JacobiBatchDerivation( degree, a, b, points );
    for ( auto const* matrix : { &values, &gradients, &runtimeValues, &runtimeGradients } )
    {
        BOOST_REQUIRE_EQUAL( matrix->size1(), degree + 1 );
        BOOST_REQUIRE_EQUAL( matrix->size2(), batchSize );
    }
    for ( unsigned q = 0; q < batchSize; ++q )
        BOOST_TEST_CONTEXT( "point=" << q << " x=" << points[q] )
        {
            checkScalar( dyna::Jacobi<double>{ static_cast<uint16_type>( degree ), a, b },
                         degree, a, b, points[q] );
            for ( unsigned i = 0; i <= degree; ++i )
            {
                auto const expected = polynomial_test::jacobi( i, a, b, points[q] );
                checkNear( values( i, q ), expected.value );
                checkNear( runtimeValues( i, q ), expected.value );
                checkNear( gradients( i, q ), expected.gradient[0] );
                checkNear( runtimeGradients( i, q ), expected.gradient[0] );
            }
        }
}

/** Compile-time Jacobi degrees dispatched by Boost.Test template cases. */
using static_orders = boost::mpl::list<boost::mpl::int_<0>, boost::mpl::int_<1>,
    boost::mpl::int_<2>, boost::mpl::int_<3>, boost::mpl::int_<5>, boost::mpl::int_<10>>;

/** @test Static Jacobi instantiations share scalar value/derivative assertions. */
BOOST_AUTO_TEST_CASE_TEMPLATE( static_jacobi, Order, static_orders )
{
    for ( std::size_t parameter = 0; parameter < alpha.size(); ++parameter )
        for ( double x : samples )
            BOOST_TEST_CONTEXT( "alpha=" << alpha[parameter] << " beta=" << beta[parameter] << " x=" << x )
                checkScalar( Jacobi<Order::value>{ alpha[parameter], beta[parameter] },
                             Order::value, alpha[parameter], beta[parameter], x );
}

/** @test Simplex values/gradients agree for each degree and normalization policy. */
BOOST_DATA_TEST_CASE( dubiner_against_polynomial_reference,
    bdata::make( dimensions ) * bdata::make( simplexDegrees ) * bdata::make( normalization ),
    dimension, degree, normalized )
{
    polynomial_test::visitDimension( dimension, [&]<int Dim>()
    {
        if ( normalized )
            checkDubiner<Dim, true>( degree );
        else
            checkDubiner<Dim, false>( degree );
    } );
}

/** @test Higher-degree triangular modes remain accurate at boundary limits. */
BOOST_DATA_TEST_CASE( high_degree_triangle_reference,
    bdata::make( std::array<unsigned, 2>{ 8, 10 } ) * bdata::make( normalization ), degree, normalized )
{
    if ( normalized )
        checkDubiner<2, true>( degree );
    else
        checkDubiner<2, false>( degree );
}

/** @test Low-degree simplex Hessians match product-rule reference derivatives. */
BOOST_DATA_TEST_CASE( simplex_hessians_against_polynomial_reference,
    bdata::make( dimensions ) * bdata::make( std::array<unsigned, 4>{ 0, 1, 2, 3 } ),
    dimension, degree )
{
    polynomial_test::visitDimension( dimension, [&]<int Dim>()
    {
        using modal_type = Feel::detail::OrthonormalPolynomialSet<
            Dim, Dynamic, Dim, Scalar, double, 0, Simplex>;
        modal_type modal{ RuntimeOrder{ static_cast<uint16_type>( degree ) } };
        typename modal_type::points_type points( Dim, Dim + 2 );
        for ( int q = 0; q < Dim + 2; ++q )
            for ( int d = 0; d < Dim; ++d )
                points( d, q ) = q == Dim + 1 ? -1.0 + 2.0 / ( Dim + 1 ) : -1.0;
        for ( int d = 0; d < Dim; ++d )
            points( d, d + 1 ) = 1.0;
        auto const hessians = modal.derivate2( points );
        BOOST_REQUIRE_EQUAL( hessians.size1(), Dim );
        BOOST_REQUIRE_EQUAL( hessians.size2(), Dim );
        for ( int q = 0; q < Dim + 2; ++q )
        {
            std::array<double, Dim> point;
            for ( int d = 0; d < Dim; ++d )
                point[d] = points( d, q );
            auto const expected = polynomial_test::dubiner<Dim>( degree, point, true );
            for ( int i = 0; i < Dim; ++i )
                for ( int j = 0; j < Dim; ++j )
                {
                    BOOST_REQUIRE_EQUAL( hessians( i, j ).size1(), expected.size() );
                    BOOST_REQUIRE_EQUAL( hessians( i, j ).size2(), points.size2() );
                    BOOST_TEST_CONTEXT( "sample=" << q << " Hessian=(" << i << "," << j << ")" )
                        for ( std::size_t mode = 0; mode < expected.size(); ++mode )
                            checkNear( hessians( i, j )( mode, q ), expected[mode].hessian[i * Dim + j], 1e-9 );
                }
        }
    } );
}

BOOST_AUTO_TEST_SUITE_END()
