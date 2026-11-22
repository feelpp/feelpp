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
#include <boost/mp11/list.hpp>

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
BOOST_DATA_TEST_CASE( batchAndRuntimeJacobi,
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
using static_orders = boost::mp11::mp_list_c<int, 0, 1, 2, 3, 5, 10>;

/** @test Static Jacobi instantiations share scalar value/derivative assertions. */
BOOST_AUTO_TEST_CASE_TEMPLATE( staticJacobi, Order, static_orders )
{
    for ( std::size_t parameter = 0; parameter < alpha.size(); ++parameter )
        for ( double x : samples )
            BOOST_TEST_CONTEXT( "alpha=" << alpha[parameter] << " beta=" << beta[parameter] << " x=" << x )
                checkScalar( Jacobi<Order::value>{ alpha[parameter], beta[parameter] },
                             Order::value, alpha[parameter], beta[parameter], x );
}

/** @test Native Eigen batches retain static degree and share the runtime recurrence. */
BOOST_AUTO_TEST_CASE_TEMPLATE( eigenJacobiStaticRuntime, Order, static_orders )
{
    constexpr int degree = Order::value;
    JacobiPolynomial<degree> fixed{ 0.25, 0.75 };
    JacobiPolynomial<Dynamic> runtime{ RuntimeOrder{ degree }, 0.25, 0.75 };
    Eigen::Matrix<double, 7, 1> points;
    for ( int q = 0; q < 7; ++q ) points( q ) = samples[q];
    auto values = fixed.values( points );
    auto derivatives = runtime.derivatives( points );
    auto seconds = runtime.derivatives<2>( points );
    static_assert( decltype( values )::rowsAtCompileTime == degree + 1 );
    static_assert( decltype( values )::colsAtCompileTime == 7 );
    static_assert( decltype( values )::usesInlineStorage );
    static_assert( decltype( derivatives )::rowsAtCompileTime == Dynamic );
    static_assert( decltype( derivatives )::colsAtCompileTime == 7 );
    static_assert( !decltype( derivatives )::usesInlineStorage );
    for ( int q = 0; q < 7; ++q )
        for ( int n = 0; n <= degree; ++n )
        {
            auto expected = polynomial_test::jacobi( n, 0.25, 0.75, points( q ) );
            checkNear( values.matrix()( n, q ), expected.value );
            checkNear( derivatives.matrix()( n, q ), expected.gradient[0] );
            checkNear( seconds.matrix()( n, q ), expected.hessian[0] );
        }
    for ( int q = 0; q < 7; ++q )
        checkNear( fixed.template derivative<2>( points( q ) ),
                   polynomial_test::jacobi( degree, 0.25, 0.75, points( q ) ).hessian[0] );
    auto const* coefficients = values.matrix().data();
    fixed.valuesInto( values, points );
    BOOST_CHECK( values.matrix().data() == coefficients );
    fixed.derivativesInto( values, points );
    BOOST_CHECK( values.matrix().data() == coefficients );
    BOOST_CHECK_SMALL( ( values.matrix() - derivatives.matrix() ).norm(), 1e-12 );
    Eigen::RowVectorXd empty( 0 );
    auto noSamples = fixed.values( empty );
    BOOST_CHECK_EQUAL( noSamples.rows(), degree + 1 );
    BOOST_CHECK_EQUAL( noSamples.cols(), 0 );
}

/** @test Strided coordinates and aliased rows are evaluated before output replacement. */
BOOST_AUTO_TEST_CASE( eigenJacobiStridesAndAliasing )
{
    JacobiPolynomial<Dynamic> polynomial{ RuntimeOrder{ 5 }, 0.25, 0.75 };
    Eigen::Matrix<double, 3, 7, Eigen::ColMajor> coordinates =
        Eigen::Matrix<double, 3, 7, Eigen::ColMajor>::Zero();
    for ( int q = 0; q < 7; ++q ) coordinates( 1, q ) = samples[q];
    auto values = polynomial.values( coordinates.row( 1 ) );
    auto expected = values.matrix().eval();
    values.matrix().row( 0 ) = coordinates.row( 1 );
    polynomial.valuesInto( values, values.matrix().row( 0 ) );
    BOOST_CHECK_SMALL( ( values.matrix() - expected ).norm(), 1e-12 );
    auto expressionValues = polynomial.values( ( 2. * coordinates.row( 1 ) ).eval() / 2. );
    BOOST_CHECK_SMALL( ( expressionValues.matrix() - expected ).norm(), 1e-12 );
    BOOST_CHECK_THROW( (void)polynomial.values( Eigen::MatrixXd::Zero( 2, 3 ) ), std::invalid_argument );
    ublas::vector<double> legacyPoints( 1 );
    BOOST_CHECK_THROW( (void)JacobiBatchEvaluation( -1, 0., 0., legacyPoints ), std::invalid_argument );
}

/** @test Mixed compile-time vector axes preserve orientation and snapshot aliased input. */
BOOST_AUTO_TEST_CASE( mixedExtentJacobiViews )
{
    JacobiPolynomial<3> polynomial;
    Eigen::RowVector3d points;
    points << -.7, -.2, .3;
    auto table = polynomial.values( points );
    auto expected = table.matrix().eval();
    table.matrix().row( 0 ) = points;
    polynomial.valuesInto( table, table.matrix().topRows( 1 ) );
    BOOST_CHECK_SMALL( ( table.matrix() - expected ).norm(), 1e-12 );
    Eigen::Matrix<double, 3, Eigen::Dynamic> columns( 3, 2 );
    columns.col( 0 ) = points.transpose();
    auto columnValues = polynomial.values( columns.leftCols( 1 ) );
    static_assert( decltype( columnValues )::colsAtCompileTime == 3 );
    BOOST_CHECK_SMALL( ( columnValues.matrix() - expected ).norm(), 1e-12 );
}

/** @test Caller-owned destinations preserve strided aliases, derivative zeros and shape checks. */
BOOST_AUTO_TEST_CASE_TEMPLATE( borrowedPolynomialDestinations, Order, static_orders )
{
    constexpr int degree = Order::value;
    Eigen::RowVector3d points;
    points << -.7, -.2, .3;
    auto checkValues = [&]( auto const& kernel )
    {
        auto expected = kernel.values( points );
        Eigen::Matrix<double, Dynamic, Dynamic, Eigen::RowMajor> storage( degree + 1, 11 );
        storage.setConstant( -99. );
        auto output = storage.block( 0, 2, degree + 1, 3 );
        int const aliasRow = degree > 1 ? degree - 1 : 0;
        output.row( aliasRow ) = points;
        kernel.valuesInto( output, output.row( aliasRow ) );
        BOOST_CHECK_SMALL( ( output - expected.matrix() ).norm(), 1e-12 );
        BOOST_CHECK_EQUAL( storage( 0, 0 ), -99. );
        BOOST_CHECK_EQUAL( storage( degree, 10 ), -99. );
        auto wrong = storage.topRows( degree );
        BOOST_CHECK_THROW( kernel.valuesInto( wrong, points ), std::invalid_argument );
    };
    JacobiPolynomial<degree> polynomial{ .25, .75 };
    checkValues( polynomial );
    checkValues( ScalingPowers<degree>{} );
    Eigen::MatrixXd derivatives = Eigen::MatrixXd::Constant( degree + 1, 3, -99. );
    polynomial.template derivativesInto<2>( derivatives, points );
    auto expected = polynomial.template derivatives<2>( points );
    BOOST_CHECK_SMALL( ( derivatives - expected.matrix() ).norm(), 1e-12 );
    BOOST_CHECK_SMALL( derivatives.topRows( std::min( degree + 1, 2 ) ).norm(), 1e-30 );
}

/** @test Static/runtime scaling tables use one recurrence for all batch sizes. */
BOOST_AUTO_TEST_CASE_TEMPLATE( scalingPowersStaticRuntime, Order, static_orders )
{
    ScalingPowers<Order::value> fixed;
    ScalingPowers<Dynamic> runtime{ RuntimeOrder{ Order::value } };
    for ( auto size : batchSizes )
    {
        Eigen::RowVectorXd points( size );
        for ( unsigned q = 0; q < size; ++q ) points( q ) = samples[q];
        auto values = fixed.values( points );
        auto dynamicValues = runtime.values( points );
        for ( int n = 0; n <= Order::value; ++n )
            for ( unsigned q = 0; q < size; ++q )
                BOOST_CHECK_SMALL( values.matrix()( n, q ) - std::pow( ( 1. - points( q ) ) / 2., n ), 1e-13 );
        BOOST_CHECK_SMALL( ( values.matrix() - dynamicValues.matrix() ).norm(), 1e-13 );
        if ( size != 0 )
        {
            auto expected = values.matrix().eval();
            values.matrix().row( 0 ) = points;
            fixed.valuesInto( values, values.matrix().row( 0 ) );
            BOOST_CHECK_SMALL( ( values.matrix() - expected ).norm(), 1e-13 );
        }
    }
}

/** Triangle/tetrahedron coordinate-map datasets. */
using coordinateDimensions = boost::mp11::mp_list_c<int, 2, 3>;

/** @test Coordinate maps preserve fixed spatial/sample extents and legacy limits. */
BOOST_AUTO_TEST_CASE_TEMPLATE( simplexCoordinateMaps, Dimension, coordinateDimensions )
{
    constexpr int dim = Dimension::value;
    Eigen::Matrix<double, dim, 3> collapsed = Eigen::Matrix<double, dim, 3>::Zero();
    for ( int d = 0; d < dim; ++d )
    {
        collapsed( d, 1 ) = -0.7 + 0.2 * d;
        collapsed( d, 2 ) = 0.1 + 0.2 * d;
    }
    auto cartesian = SimplexCoordinates<dim>::toCartesian( collapsed );
    auto restored = SimplexCoordinates<dim>::toCollapsed( cartesian.matrix() );
    static_assert( decltype( cartesian )::rowsAtCompileTime == dim );
    static_assert( decltype( cartesian )::colsAtCompileTime == 3 );
    BOOST_CHECK_SMALL( ( restored.matrix() - collapsed ).norm(), 1e-13 );
    BOOST_CHECK_EQUAL( cartesian.matrix()( 0, 0 ), dim == 2 ? -.5 : -.75 );
    BOOST_CHECK_EQUAL( cartesian.matrix()( dim - 1, 0 ), 0. );
    Eigen::Matrix<double, dim, 1> vertex = Eigen::Matrix<double, dim, 1>::Constant( -1. );
    vertex( dim - 1 ) = 1.;
    auto representative = SimplexCoordinates<dim>::toCollapsed( vertex );
    static_assert( decltype( representative )::usesInlineStorage );
    BOOST_CHECK_EQUAL( representative.matrix()( 0, 0 ), dim == 2 ? -1. : 1. );
    if constexpr ( dim == 3 ) BOOST_CHECK_EQUAL( representative.matrix()( 1, 0 ), -1. );
    BOOST_CHECK_EQUAL( representative.matrix()( dim - 1, 0 ), 1. );
    constexpr TheShape shape = dim == 2 ? TRIANGLE : TETRAHEDRON;
    typename node<double>::type legacy( dim );
    eigenMap( legacy ) = vertex;
    details::CollapsedCoordinates<shape> map;
    details::CollapsedCoordinates<shape> constructed( legacy );
    BOOST_CHECK_SMALL( ( eigenMap( map( legacy ) ) - representative.matrix().col( 0 ) ).norm(), 1e-13 );
    BOOST_CHECK_SMALL( ( eigenMap( constructed() ) - eigenMap( map() ) ).norm(), 1e-13 );
    for ( int count : { 0, 1, 3 } )
    {
        typename matrix_node<double>::type batch( dim, count );
        eigenMap( batch ) = cartesian.matrix().leftCols( count );
        details::CollapsedPointCoordinates<shape> batchMap;
        auto expected = SimplexCoordinates<dim>::toCollapsed( eigenMap( batch ) );
        BOOST_CHECK_SMALL( ( eigenMap( batchMap( batch ) ) - expected.matrix() ).norm(), 1e-13 );
        auto aliasedExpected = SimplexCoordinates<dim>::toCollapsed( eigenMap( batchMap() ) );
        BOOST_CHECK_SMALL( ( eigenMap( batchMap( batchMap() ) ) - aliasedExpected.matrix() ).norm(), 1e-13 );
    }
}

/** @test Simplex values/gradients agree for each degree and normalization policy. */
BOOST_DATA_TEST_CASE( dubinerAgainstPolynomialReference,
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
BOOST_DATA_TEST_CASE( highDegreeTriangleReference,
    bdata::make( std::array<unsigned, 2>{ 8, 10 } ) * bdata::make( normalization ), degree, normalized )
{
    if ( normalized )
        checkDubiner<2, true>( degree );
    else
        checkDubiner<2, false>( degree );
}

/** @test Low-degree simplex Hessians match product-rule reference derivatives. */
BOOST_DATA_TEST_CASE( simplexHessiansAgainstPolynomialReference,
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
