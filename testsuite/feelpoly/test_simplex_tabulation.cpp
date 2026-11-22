/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file test_simplex_tabulation.cpp
 * @brief Independent boundary jets, ordering, orthogonality and alias contracts.
 */
#define BOOST_TEST_MODULE simplex_tabulation
#include <feel/feelcore/testsuite.hpp>
#include <feel/feelpoly/simplextabulation.hpp>
#include <feel/feelpoly/dubiner.hpp>
#include <boost/mp11/list.hpp>
#include <Eigen/Eigenvalues>
#include "polynomial_reference.hpp"
#include <array>
#include <cmath>
#include <exception>
#include <thread>
#include <vector>

using namespace Feel;
FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Dimensions exercise the same reference identities and ownership contract. */
using Dimensions = boost::mp11::mp_list_c<int, 1, 2, 3>;

/** @brief Compare an actual derivative with the independent 100-digit oracle. */
void checkReference( double actual, polynomial_test::precise_type const& expected )
{
    double const reference = expected.convert_to<double>();
    BOOST_REQUIRE( std::isfinite( actual ) );
    BOOST_CHECK_SMALL( actual - reference, 3e-11 * ( 1. + std::abs( reference ) ) );
}

/** @brief Sample every entity and two interior approaches to every vertex.
 * Entity centers are generated independently from barycentric vertex subsets.
 */
template<int Dim>
auto boundaryPoints()
{
    using Point = Eigen::Matrix<double, Dim, 1>;
    std::array<Point, Dim + 1> vertices;
    for ( auto& vertex : vertices ) vertex.setConstant( -1. );
    for ( int d = 0; d < Dim; ++d ) vertices[d + 1]( d ) = 1.;
    std::vector<Point> samples;
    for ( unsigned mask = 1; mask < ( 1u << ( Dim + 1 ) ); ++mask )
    {
        Point point = Point::Zero();
        int count = 0;
        for ( int v = 0; v <= Dim; ++v )
            if ( mask & ( 1u << v ) ) { point += vertices[v]; ++count; }
        samples.push_back( point / count );
    }
    std::array<Point, 2> interiors;
    interiors[0].setConstant( -1. + 2. / ( Dim + 1 ) );
    interiors[1].setZero();
    for ( int v = 0; v <= Dim; ++v )
        interiors[1] += ( 2. * ( v + 1 ) / ( ( Dim + 1 ) * ( Dim + 2 ) ) ) * vertices[v];
    for ( auto const& vertex : vertices )
        for ( auto const& interior : interiors )
            for ( double epsilon : { 1e-4, 1e-9, 1e-14 } )
                samples.push_back( ( 1. - epsilon ) * vertex + epsilon * interior );
    Eigen::Matrix<double, Dim, Dynamic> points( Dim, samples.size() );
    for ( std::size_t q = 0; q < samples.size(); ++q ) points.col( q ) = samples[q];
    return points;
}

/** @brief Check both order policies against finite binomial expansion and product rules. */
template<int Dim, int P, bool Normalized>
void checkBoundaryJets()
{
    using Static = SimplexDubinerBasis<Dim, P, double, Normalized>;
    using Runtime = SimplexDubinerBasis<Dim, Dynamic, double, Normalized>;
    Static fixed;
    Runtime runtime{ RuntimeOrder{ P } };
    auto points = boundaryPoints<Dim>();
    auto values = fixed.values( points );
    auto gradients = fixed.gradients( points );
    auto hessians = fixed.hessians( points );
    auto dynamicValues = runtime.values( points );
    auto dynamicGradients = runtime.gradients( points );
    auto dynamicHessians = runtime.hessians( points );
    static_assert( decltype( values )::rowsAtCompileTime == Static::basisExtent );
    static_assert( decltype( dynamicValues )::rowsAtCompileTime == Dynamic );
    static_assert( decltype( values )::colsAtCompileTime == Dynamic );
    BOOST_CHECK_SMALL( ( values.matrix() - dynamicValues.matrix() ).norm(), 1e-11 );
    for ( int d = 0; d < Dim; ++d )
        BOOST_CHECK_SMALL( ( gradients[d].matrix() - dynamicGradients[d].matrix() ).norm(), 1e-10 );
    for ( int ij = 0; ij < Dim * Dim; ++ij )
        BOOST_CHECK_SMALL( ( hessians[ij].matrix() - dynamicHessians[ij].matrix() ).norm(), 1e-9 );
    for ( Eigen::Index q = 0; q < points.cols(); ++q )
    {
        std::array<double, Dim> point;
        for ( int d = 0; d < Dim; ++d ) point[d] = points( d, q );
        auto const expected = polynomial_test::dubiner<Dim>( P, point, Normalized );
        BOOST_REQUIRE_EQUAL( expected.size(), values.rows() );
        BOOST_TEST_CONTEXT( "dimension=" << Dim << " degree=" << P << " normalized=" << Normalized << " sample=" << q )
            for ( std::size_t row = 0; row < expected.size(); ++row )
            {
                checkReference( values.matrix()( row, q ), expected[row].value );
                for ( int i = 0; i < Dim; ++i )
                {
                    checkReference( gradients[i].matrix()( row, q ), expected[row].gradient[i] );
                    for ( int j = 0; j < Dim; ++j )
                    {
                        checkReference( hessians[i * Dim + j].matrix()( row, q ), expected[row].hessian[i * Dim + j] );
                        BOOST_CHECK_EQUAL( hessians[i * Dim + j].matrix()( row, q ), hessians[j * Dim + i].matrix()( row, q ) );
                    }
                }
            }
    }
    Eigen::Index row = 0;
    for ( int total = 0; total <= P; ++total )
        for ( int k = 0; k <= ( Dim == 3 ? total : 0 ); ++k )
            for ( int j = 0; j <= ( Dim >= 2 ? total - k : 0 ); ++j )
            {
                typename Static::MultiIndex degrees{};
                degrees[0] = total - j - k;
                if constexpr ( Dim >= 2 ) degrees[1] = j;
                if constexpr ( Dim == 3 ) degrees[2] = k;
                BOOST_CHECK_EQUAL( fixed.modeIndex( degrees ), row++ );
            }
    BOOST_CHECK_EQUAL( row, fixed.basisCount() );
    typename Static::MultiIndex invalid{};
    invalid[0] = P + 1;
    BOOST_CHECK_THROW( (void)fixed.modeIndex( invalid ), std::out_of_range );
    invalid[0] = -1;
    BOOST_CHECK_THROW( (void)fixed.modeIndex( invalid ), std::out_of_range );
    Eigen::Matrix<double, Dim, Dynamic> empty( Dim, 0 );
    BOOST_CHECK_EQUAL( fixed.values( empty ).cols(), 0 );
    BOOST_CHECK_EQUAL( runtime.gradients( empty )[0].cols(), 0 );
    BOOST_CHECK_EQUAL( runtime.hessians( empty )[0].cols(), 0 );
    BOOST_CHECK_THROW( (void)fixed.values( Eigen::MatrixXd::Zero( Dim + 1, 2 ) ), std::invalid_argument );
}

/** @brief Gauss-Legendre nodes/weights from an independent symmetric eigenproblem. */
auto gaussLegendre( int count )
{
    Eigen::MatrixXd jacobi = Eigen::MatrixXd::Zero( count, count );
    for ( int n = 1; n < count; ++n )
        jacobi( n - 1, n ) = jacobi( n, n - 1 ) = n / std::sqrt( 4. * n * n - 1. );
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver( jacobi );
    if ( solver.info() != Eigen::Success ) throw std::runtime_error( "test quadrature eigensolve failed" );
    Eigen::VectorXd nodes = solver.eigenvalues();
    Eigen::VectorXd weights = 2. * solver.eigenvectors().row( 0 ).array().square().transpose();
    return std::pair{ nodes, weights };
}

/** @brief Integrate the complete modal Gram matrix with exact polynomial cubature. */
template<int Dim>
void checkGramMatrix( int degree )
{
    auto const [nodes, weights] = gaussLegendre( degree + 2 );
    int const count = degree + 2;
    int samples = 1;
    for ( int d = 0; d < Dim; ++d ) samples *= count;
    Eigen::Matrix<double, Dim, Dynamic> points( Dim, samples );
    Eigen::VectorXd quadrature( samples );
    for ( int q = 0; q < samples; ++q )
    {
        int index = q;
        std::array<double, Dim> coordinates;
        double weight = 1.;
        for ( int d = 0; d < Dim; ++d )
        {
            coordinates[d] = nodes( index % count );
            weight *= weights( index % count );
            index /= count;
        }
        double factor = 1.;
        for ( int d = Dim - 1; d >= 0; --d )
        {
            points( d, q ) = ( 1. + coordinates[d] ) * factor - 1.;
            weight *= factor;
            factor *= ( 1. - coordinates[d] ) / 2.;
        }
        quadrature( q ) = weight;
    }
    SimplexDubinerBasis<Dim, Dynamic> basis{ RuntimeOrder{ degree } };
    auto values = basis.values( points );
    Eigen::MatrixXd gram = values.matrix() * quadrature.asDiagonal() * values.matrix().transpose();
    Eigen::MatrixXd expected = Eigen::MatrixXd::Identity( gram.rows(), gram.cols() );
    BOOST_CHECK_SMALL( ( gram - expected ).cwiseAbs().maxCoeff(), 3e-12 );
    double volume = 2.;
    if constexpr ( Dim == 3 ) volume = 4. / 3.;
    BOOST_CHECK_SMALL( quadrature.sum() - volume, 1e-13 );
    BOOST_CHECK_SMALL( values.matrix()( 0, 0 ) - 1. / std::sqrt( volume ), 1e-14 );
}
} // namespace

/** @test All boundary strata and vertex approaches use the same regular polynomial path. */
BOOST_AUTO_TEST_CASE_TEMPLATE( simplexBoundaryJets, Dimension, Dimensions )
{
    constexpr int dim = Dimension::value;
    checkBoundaryJets<dim, 0, false>();
    checkBoundaryJets<dim, 1, false>();
    checkBoundaryJets<dim, 2, false>();
    checkBoundaryJets<dim, 3, false>();
    checkBoundaryJets<dim, 5, false>();
    checkBoundaryJets<dim, 0, true>();
    checkBoundaryJets<dim, 1, true>();
    checkBoundaryJets<dim, 2, true>();
    checkBoundaryJets<dim, 3, true>();
    checkBoundaryJets<dim, 5, true>();
}

/** @test Orthonormality, constant normalization and complete cardinality agree in all dimensions. */
BOOST_AUTO_TEST_CASE_TEMPLATE( simplexGramMatrix, Dimension, Dimensions )
{
    for ( int degree : { 0, 1, 2, 3, 5 } )
        BOOST_TEST_CONTEXT( "degree=" << degree ) { checkGramMatrix<Dimension::value>( degree ); }
}

/** @test Reused storage handles strided/expression points and aliases in every derivative result. */
BOOST_AUTO_TEST_CASE_TEMPLATE( simplexWorkspaceOwnership, Dimension, Dimensions )
{
    constexpr int dim = Dimension::value;
    using Basis = SimplexDubinerBasis<dim, 3>;
    Basis basis;
    Eigen::Matrix<double, dim, 3> points;
    points.col( 0 ).setConstant( -1. );
    points.col( 1 ).setConstant( -1. + 2. / ( dim + 1 ) );
    points.col( 2 ).setConstant( -1. );
    points( dim - 1, 2 ) = 1.;
    auto values = basis.values( points );
    auto gradients = basis.gradients( points );
    auto hessians = basis.hessians( points );
    static_assert( decltype( values )::colsAtCompileTime == 3 );
    static_assert( decltype( values )::usesInlineStorage );
    auto checkAlias = [&]( auto& output, auto& coordinateOwner, auto& workspace, auto&& evaluate )
    {
        auto expected = output;
        auto const* address = coordinateOwner.matrix().data();
        coordinateOwner.matrix().topRows( dim ) = points;
        evaluate( output, coordinateOwner.matrix().topRows( dim ), workspace );
        BOOST_CHECK_EQUAL( coordinateOwner.matrix().data(), address );
        if constexpr ( requires { output.matrix(); } )
            BOOST_CHECK_SMALL( ( output.matrix() - expected.matrix() ).norm(), 1e-12 );
        else
            for ( std::size_t d = 0; d < output.size(); ++d )
                BOOST_CHECK_SMALL( ( output[d].matrix() - expected[d].matrix() ).norm(), 1e-11 );
    };
    typename Basis::template Workspace<0> valueWorkspace;
    typename Basis::template Workspace<1> gradientWorkspace;
    typename Basis::template Workspace<2> hessianWorkspace;
    checkAlias( values, values, valueWorkspace, [&]( auto& output, auto const& input, auto& workspace )
    { basis.valuesInto( output, input, workspace ); } );
    checkAlias( gradients, gradients.back(), gradientWorkspace, [&]( auto& output, auto const& input, auto& workspace )
    { basis.gradientsInto( output, input, workspace ); } );
    checkAlias( hessians, hessians.back(), hessianWorkspace, [&]( auto& output, auto const& input, auto& workspace )
    { basis.hessiansInto( output, input, workspace ); } );
    auto expression = basis.hessians( points + points * 0. );
    Eigen::Matrix<double, dim + 1, 3> storage;
    storage.setZero();
    storage.bottomRows( dim ) = points;
    auto strided = basis.hessians( storage.bottomRows( dim ) );
    for ( int ij = 0; ij < dim * dim; ++ij )
    {
        BOOST_CHECK_SMALL( ( expression[ij].matrix() - hessians[ij].matrix() ).norm(), 1e-12 );
        BOOST_CHECK_SMALL( ( strided[ij].matrix() - hessians[ij].matrix() ).norm(), 1e-12 );
    }
    SimplexDubinerBasis<dim, Dynamic> runtime{ RuntimeOrder{ 3 } };
    typename decltype( runtime )::template Workspace<0> wrong{ RuntimeOrder{ 2 } };
    auto output = runtime.makeValues( points.cols() );
    BOOST_CHECK_THROW( runtime.valuesInto( output, points, wrong ), std::invalid_argument );
}

/** @test Preserve tiny affine numerators before high-degree derivative amplification.
 * This P16 boundary witness is checked against finite binomial expansion in
 * 100 digits, independently of either production or external double recurrences.
 */
BOOST_AUTO_TEST_CASE( smallBoundaryDerivativeRetainsPrecision )
{
    Eigen::Vector3d point;
    point << -.9999999995000001, -.9999999995000001, .9999999985;
    SimplexDubinerBasis<3, Dynamic> basis{ RuntimeOrder{ 16 } };
    auto hessians = basis.hessians( point );
    using Precise = polynomial_test::precise_type;
    Precise const u = Precise( 1 ) + Precise( point( 0 ) ) +
        ( Precise( point( 1 ) ) + Precise( point( 2 ) ) ) / 2;
    auto const jacobi = polynomial_test::jacobi( 13, 8., 0., point( 2 ) );
    Precise const expected = 15 * u * jacobi.value * sqrt( Precise( 245 ) );
    auto const row = basis.modeIndex( { 3, 0, 13 } );
    BOOST_CHECK_SMALL( hessians[0].matrix()( row, 0 ) - expected.convert_to<double>(), 1e-20 );
}

/** @test Legacy values and coefficient operators agree with direct sampled derivatives. */
BOOST_AUTO_TEST_CASE_TEMPLATE( dubinerCompatibilityAndOperatorIdentities, Dimension, Dimensions )
{
    constexpr int dim = Dimension::value;
    using Basis = Dubiner<dim, dim, 3>;
    Basis basis;
    auto points = boundaryPoints<dim>();
    typename Basis::points_type legacy( dim, points.cols() );
    eigenMap( legacy ) = points;
    auto values = Basis::values( points );
    auto gradients = Basis::gradients( points );
    auto hessians = Basis::hessians( points );
    auto compatibility = Basis::evaluate( legacy );
    BOOST_CHECK_SMALL( ( eigenMap( compatibility ) - values.matrix() ).norm(), 1e-12 );
    auto copy = basis;
    for ( int i = 0; i < dim; ++i )
    {
        BOOST_CHECK_EQUAL( &copy.d( i ), &basis.d( i ) );
        auto differentiated = ( eigenMap( basis.d( i ) ) * values.matrix() ).eval();
        BOOST_CHECK_SMALL( ( differentiated - gradients[i].matrix() ).cwiseAbs().maxCoeff(), 1e-10 );
        for ( int j = 0; j < dim; ++j )
            BOOST_CHECK_SMALL( ( eigenMap( basis.d( i ) ) * eigenMap( basis.d( j ) ) * values.matrix() -
                hessians[i * dim + j].matrix() ).cwiseAbs().maxCoeff(), 1e-9 );
    }
    BOOST_CHECK_THROW( (void)basis.d( dim ), std::out_of_range );
    using Native = Dubiner<dim, dim, 3, Normalized<true>, double, StorageEigen>;
    Native native;
    auto nativeValues = Native::evaluate( points );
    BOOST_CHECK_SMALL( ( nativeValues - values.matrix() ).norm(), 1e-12 );
    BOOST_CHECK_SMALL( ( native.d( 0 ) * nativeValues - gradients[0].matrix() ).cwiseAbs().maxCoeff(), 1e-10 );
}

/** @test Concurrent first operator requests share one immutable, lifetime-stable result. */
BOOST_AUTO_TEST_CASE( concurrentDubinerOperatorInitialization )
{
    using Basis = Dubiner<2, 2, 4, Normalized<false>>;
    std::array<typename Basis::matrix_type const*, 4> addresses{};
    std::array<std::exception_ptr, 4> errors{};
    std::array<std::thread, 4> workers;
    for ( std::size_t i = 0; i < workers.size(); ++i )
        workers[i] = std::thread( [&, i]
        {
            try { Basis basis; addresses[i] = &basis.d( 0 ); }
            catch ( ... ) { errors[i] = std::current_exception(); }
        } );
    for ( auto& worker : workers ) worker.join();
    for ( std::size_t i = 0; i < workers.size(); ++i )
    {
        BOOST_REQUIRE( !errors[i] );
        BOOST_REQUIRE( addresses[i] );
        BOOST_CHECK_EQUAL( addresses[i], addresses[0] );
    }
    BOOST_CHECK_EQUAL( addresses[0]->size1(), 15 );
}

/** @test Value-only/high-order construction does not require an auxiliary operator solve. */
BOOST_AUTO_TEST_CASE( highOrderDubinerValuesWithoutOperatorSetup )
{
    constexpr Dubiner<2, 2, 20> basis;
    Eigen::Matrix<double, 2, 1> vertex;
    vertex << -1., 1.;
    auto values = basis.values( vertex );
    static_assert( decltype( values )::rowsAtCompileTime == 231 );
    BOOST_CHECK_EQUAL( values.rows(), 231 );
    BOOST_CHECK( values.matrix().allFinite() );
    auto hessians = basis.hessians( vertex );
    for ( auto const& entry : hessians ) BOOST_CHECK( entry.matrix().allFinite() );
}
