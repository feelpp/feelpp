/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_geomap_inverse.cpp
 * @brief Structured geometric inverse and embedded closest-point tests.
 */
#define BOOST_TEST_MODULE test_geomap_inverse
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/geomap.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Geometric coefficient owner satisfying the `GeoMap::Inverse` constructor. */
template <typename GeomapType>
struct synthetic_geometry
{
    using matrix_type = typename GeomapType::matrix_node_t_type;

    /** Physical interpolation coefficients. */
    matrix_type coefficients;

    /** @brief Return all geometric coefficients. */
    [[nodiscard]] matrix_type const& G() const noexcept { return coefficients; }
    /** @brief Return coefficients for the affine compatibility constructor. */
    [[nodiscard]] matrix_type const& vertices() const noexcept
    {
        return coefficients;
    }
};

/** Build interpolation coefficients by sampling a polynomial physical map. */
template <typename GeomapType, typename Mapping>
[[nodiscard]] synthetic_geometry<GeomapType>
makeGeometry( std::shared_ptr<GeomapType> const& geomap, Mapping&& mapping )
{
    synthetic_geometry<GeomapType> geometry{
        typename GeomapType::matrix_node_t_type(
            GeomapType::nRealDim, geomap->nbPoints() ) };
    for ( uint16_type point = 0; point < geomap->nbPoints(); ++point )
    {
        Eigen::Matrix<double, GeomapType::nDim, 1> referencePoint;
        for ( int component = 0; component < GeomapType::nDim; ++component )
            referencePoint( component ) =
                geomap->refNode( point )( component );
        auto const physicalPoint = std::invoke( mapping, referencePoint );
        for ( int component = 0; component < GeomapType::nRealDim;
              ++component )
            geometry.coefficients( component, point ) =
                physicalPoint( component );
    }
    return geometry;
}

/** Convert a fixed Eigen vector to a Feel++ node. */
template <typename EigenVector>
[[nodiscard]] node_t<double> makeNode( EigenVector const& point )
{
    node_t<double> result( point.rows() );
    for ( Eigen::Index component = 0; component < point.rows(); ++component )
        result( component ) = point( component );
    return result;
}

/** Return the polynomial test mapping represented exactly at `Order`. */
template <int Order>
[[nodiscard]] Eigen::Vector2d
polynomialMap( Eigen::Vector2d const& point )
{
    double const x = point( 0 );
    double const y = point( 1 );
    Eigen::Vector2d result;
    result << 0.2 + 1.15 * x + 0.12 * y,
        -0.3 + 0.08 * x + 1.1 * y;
    if constexpr ( Order >= 2 )
    {
        result( 0 ) += 0.08 * x * y + 0.035 * y * y;
        result( 1 ) += 0.06 * x * x - 0.025 * x * y;
    }
    if constexpr ( Order >= 3 )
    {
        result( 0 ) += 0.018 * x * x * x - 0.012 * x * y * y;
        result( 1 ) += 0.015 * y * y * y + 0.01 * x * x * y;
    }
    return result;
}

/** Check exact interior and boundary round trips for one geometric order. */
template <int Order>
void checkPolynomialRoundTrips()
{
    using geomap_type = GeoMap<2, Order, 2, double, Simplex>;
    auto geomap = std::make_shared<geomap_type>();
    auto geometry = makeGeometry(
        geomap, []( Eigen::Vector2d const& point )
        { return polynomialMap<Order>( point ); } );
    typename geomap_type::Inverse inverse( geomap, geometry );

    for ( Eigen::Vector2d const referencePoint : {
              Eigen::Vector2d( -0.7, -0.55 ),
              Eigen::Vector2d( 0.15, -0.45 ),
              Eigen::Vector2d( -0.45, 0.12 ) } )
    {
        auto const result = inverse.inverse(
            makeNode( polynomialMap<Order>( referencePoint ) ) );
        BOOST_CHECK( result.converged() );
        BOOST_CHECK( result.exact() );
        BOOST_CHECK( result.location ==
                     Feel::geometry::InverseLocation::inside );
        BOOST_CHECK_SMALL(
            ( result.reference_point - referencePoint ).norm(), 2e-10 );
        BOOST_CHECK_SMALL( result.residual_norm, 2e-11 );
        BOOST_CHECK_EQUAL( result.jacobian_rank, 2 );
    }

    Eigen::Vector2d const boundaryPoint( -0.2, 0.2 );
    auto const boundaryResult = inverse.inverse(
        makeNode( polynomialMap<Order>( boundaryPoint ) ) );
    BOOST_CHECK( boundaryResult.exact() );
    BOOST_CHECK( boundaryResult.location ==
                 Feel::geometry::InverseLocation::boundary );
    BOOST_CHECK_SMALL(
        ( boundaryResult.reference_point - boundaryPoint ).norm(), 2e-10 );
}
} // namespace

BOOST_AUTO_TEST_SUITE( geomap_inverse_suite )

/** @test Verify exact full-dimensional G1, G2, and G3 round trips. */
BOOST_AUTO_TEST_CASE( polynomial_round_trips_g1_g2_g3 )
{
    checkPolynomialRoundTrips<1>();
    checkPolynomialRoundTrips<2>();
    checkPolynomialRoundTrips<3>();
}

/** @test Verify static/runtime order parity through the shared inverse kernel. */
BOOST_AUTO_TEST_CASE( static_runtime_order_parity )
{
    using static_geomap_type = GeoMap<2, 2, 2, double, Simplex>;
    using runtime_geomap_type = GeoMap<2, Dynamic, 2, double, Simplex>;
    auto staticGeomap = std::make_shared<static_geomap_type>();
    auto runtimeGeomap = std::make_shared<runtime_geomap_type>(
        RuntimeOrder{ 2 } );
    auto mapping = []( Eigen::Vector2d const& point )
    { return polynomialMap<2>( point ); };
    auto staticGeometry = makeGeometry( staticGeomap, mapping );
    auto runtimeGeometry = makeGeometry( runtimeGeomap, mapping );
    typename static_geomap_type::Inverse staticInverse(
        staticGeomap, staticGeometry );
    typename runtime_geomap_type::Inverse runtimeInverse(
        runtimeGeomap, runtimeGeometry );

    Eigen::Vector2d const referencePoint( -0.18, -0.37 );
    auto const target = makeNode( mapping( referencePoint ) );
    auto const staticResult = staticInverse.inverse( target );
    auto const runtimeResult = runtimeInverse.inverse( target );
    BOOST_REQUIRE( staticResult.converged() );
    BOOST_REQUIRE_MESSAGE(
        runtimeResult.converged(),
        "runtime inverse status=" << Feel::geometry::inverseStatusName(
                                         runtimeResult.status )
                                  << ", residual=" << runtimeResult.residual_norm
                                  << ", gradient=" << runtimeResult.gradient_norm
                                  << ", rank=" << runtimeResult.jacobian_rank
                                  << ", runtime order=" << runtimeGeomap->order() );
    BOOST_CHECK_SMALL(
        ( staticResult.reference_point -
          runtimeResult.reference_point )
            .norm(),
        2e-12 );
    BOOST_CHECK_SMALL(
        std::abs( staticResult.residual_norm - runtimeResult.residual_norm ),
        2e-13 );
}

/** @test Verify the vectorized runtime G1 batch follows static G1 exactly. */
BOOST_AUTO_TEST_CASE( affine_static_runtime_batch_parity )
{
    using static_geomap_type = GeoMap<2, 1, 2, double, Simplex>;
    using runtime_geomap_type = GeoMap<2, Dynamic, 2, double, Simplex>;
    auto staticGeomap = std::make_shared<static_geomap_type>();
    auto runtimeGeomap = std::make_shared<runtime_geomap_type>(
        RuntimeOrder{ 1 } );
    auto mapping = []( Eigen::Vector2d const& point )
    {
        Eigen::Vector2d physical;
        physical << 0.2 + 1.1 * point( 0 ) + 0.15 * point( 1 ),
            -0.3 + 0.08 * point( 0 ) + 0.95 * point( 1 );
        return physical;
    };
    auto staticGeometry = makeGeometry( staticGeomap, mapping );
    auto runtimeGeometry = makeGeometry( runtimeGeomap, mapping );
    typename static_geomap_type::Inverse staticInverse(
        staticGeomap, staticGeometry );
    typename runtime_geomap_type::Inverse runtimeInverse(
        runtimeGeomap, runtimeGeometry );

    Eigen::Matrix<double, 2, Eigen::Dynamic> targets( 2, 64 );
    std::vector<Eigen::Vector2d> referencePoints;
    referencePoints.reserve( targets.cols() );
    for ( Eigen::Index point = 0; point < targets.cols(); ++point )
    {
        double const parameter =
            double( point ) / double( targets.cols() - 1 );
        Eigen::Vector2d referencePoint(
            -0.75 + 0.35 * parameter,
            -0.7 + 0.2 * parameter );
        referencePoints.push_back( referencePoint );
        targets.col( point ) = mapping( referencePoint );
    }

    runtimeInverse.setGeometryEvaluationStrategy(
        geometry::geometry_evaluation_strategy::pointwise );
    auto const pointwiseResults = runtimeInverse.inverse( targets );
    runtimeInverse.setGeometryEvaluationStrategy(
        geometry::geometry_evaluation_strategy::batched );
    auto const staticResults = staticInverse.inverse( targets );
    auto const runtimeResults = runtimeInverse.inverse( targets );
    BOOST_REQUIRE_EQUAL( staticResults.size(), runtimeResults.size() );
    BOOST_REQUIRE_EQUAL( pointwiseResults.size(), runtimeResults.size() );
    for ( std::size_t point = 0; point < runtimeResults.size(); ++point )
    {
        BOOST_REQUIRE( staticResults[point].exact() );
        BOOST_REQUIRE( pointwiseResults[point].exact() );
        BOOST_REQUIRE( runtimeResults[point].exact() );
        BOOST_CHECK_EQUAL( staticResults[point].iterations, 0 );
        BOOST_CHECK_EQUAL( runtimeResults[point].iterations, 0 );
        BOOST_CHECK_SMALL(
            ( staticResults[point].reference_point -
              runtimeResults[point].reference_point )
                .norm(),
            2e-13 );
        BOOST_CHECK_SMALL(
            ( pointwiseResults[point].reference_point -
              runtimeResults[point].reference_point )
                .norm(),
            2e-13 );
        BOOST_CHECK_SMALL(
            ( runtimeResults[point].reference_point -
              referencePoints[point] )
                .norm(),
            2e-12 );
    }
}

/** @test Verify zero-copy views and Eigen/uBLAS batch API parity. */
BOOST_AUTO_TEST_CASE( eigen_view_and_batch_interoperability )
{
    using geomap_type = GeoMap<2, 2, 2, double, Simplex>;
    auto geomap = std::make_shared<geomap_type>();
    auto mapping = []( Eigen::Vector2d const& point )
    { return polynomialMap<2>( point ); };
    auto geometry = makeGeometry( geomap, mapping );
    typename geomap_type::Inverse inverse( geomap, geometry );

    auto coefficients = inverse.geometryCoefficientsView();
    auto jacobian = inverse.jacobianView();
    auto inverseJacobian = inverse.inverseJacobianTransposeView();
    static_assert(
        std::same_as<decltype( coefficients.data() ), double const*> );
    static_assert( std::same_as<decltype( jacobian.data() ), double const*> );
    static_assert(
        std::same_as<decltype( inverseJacobian.data() ), double const*> );
    BOOST_CHECK_EQUAL(
        coefficients.data(),
        std::to_address( inverse.G().data().begin() ) );
    BOOST_CHECK_EQUAL( jacobian.data(),
                       std::to_address( inverse.K().data().begin() ) );
    BOOST_CHECK_EQUAL(
        inverseJacobian.data(),
        std::to_address( inverse.B().data().begin() ) );

    typename geomap_type::Inverse::eigen_physical_points_type eigenTargets(
        2, 3 );
    std::array<Eigen::Vector2d, 3> const referencePoints{
        Eigen::Vector2d( -0.65, -0.55 ),
        Eigen::Vector2d( -0.45, -0.48 ),
        Eigen::Vector2d( -0.3, -0.42 ) };
    for ( std::size_t point = 0; point < referencePoints.size(); ++point )
        eigenTargets.col( point ) = mapping( referencePoints[point] );

    auto const eigenResults = inverse.inverse( eigenTargets );
    typename geomap_type::matrix_node_t_type ublasTargets( 2, 3 );
    eigenMap( ublasTargets ) = eigenTargets;
    auto const ublasResults = inverse.inverse( ublasTargets );
    BOOST_REQUIRE_EQUAL( eigenResults.size(), ublasResults.size() );
    for ( std::size_t point = 0; point < eigenResults.size(); ++point )
    {
        BOOST_REQUIRE( eigenResults[point].converged() );
        BOOST_REQUIRE( ublasResults[point].converged() );
        BOOST_CHECK_SMALL(
            ( eigenResults[point].reference_point -
              ublasResults[point].reference_point )
                .norm(),
            2e-13 );
        BOOST_CHECK_SMALL(
            ( eigenResults[point].reference_point - referencePoints[point] )
                .norm(),
            2e-10 );
    }
}

/** @test Verify unrestricted and constrained reference-domain policies. */
BOOST_AUTO_TEST_CASE( outside_and_boundary_policies )
{
    using simplex_geomap_type = GeoMap<2, 2, 2, double, Simplex>;
    auto simplexGeomap = std::make_shared<simplex_geomap_type>();
    auto simplexGeometry = makeGeometry(
        simplexGeomap, []( Eigen::Vector2d const& point )
        { return polynomialMap<2>( point ); } );
    typename simplex_geomap_type::Inverse simplexInverse(
        simplexGeomap, simplexGeometry );

    Eigen::Vector2d const outsideReference( 0.35, 0.25 );
    auto const outsideTarget = makeNode(
        polynomialMap<2>( outsideReference ) );
    auto const unrestricted = simplexInverse.inverse( outsideTarget );
    BOOST_CHECK( unrestricted.exact() );
    BOOST_CHECK( unrestricted.location ==
                 Feel::geometry::InverseLocation::outside );
    BOOST_CHECK( !simplexInverse.isIn() );

    typename simplex_geomap_type::Inverse::inverse_options_type constrained;
    constrained.domain_policy =
        Feel::geometry::InverseDomainPolicy::constrained;
    constrained.warm_start = false;
    auto const projected = simplexInverse.inverse(
        outsideTarget, constrained );
    BOOST_CHECK_MESSAGE(
        projected.converged(),
        "constrained simplex inverse status=" << Feel::geometry::inverseStatusName( projected.status )
                                              << ", residual=" << projected.residual_norm
                                              << ", gradient=" << projected.gradient_norm
                                              << ", iterations=" << projected.iterations
                                              << ", reference=" << projected.reference_point.transpose() );
    BOOST_CHECK( !projected.exact() );
    BOOST_CHECK( projected.location ==
                 Feel::geometry::InverseLocation::boundary );
    BOOST_CHECK_GT( projected.residual_norm, 1e-3 );

    using hypercube_geomap_type = GeoMap<2, 1, 2, double, Hypercube>;
    auto hypercubeGeomap = std::make_shared<hypercube_geomap_type>();
    auto hypercubeGeometry = makeGeometry(
        hypercubeGeomap, []( Eigen::Vector2d const& point )
        { return point; } );
    typename hypercube_geomap_type::Inverse hypercubeInverse(
        hypercubeGeomap, hypercubeGeometry );
    auto const hypercubeProjected = hypercubeInverse.inverse(
        makeNode( Eigen::Vector2d( 1.6, 0.3 ) ), constrained );
    BOOST_CHECK( hypercubeProjected.converged() );
    BOOST_CHECK( hypercubeProjected.location ==
                 Feel::geometry::InverseLocation::boundary );
    BOOST_CHECK_SMALL( hypercubeProjected.reference_point( 0 ) - 1.0, 1e-13 );
    BOOST_CHECK_SMALL( hypercubeProjected.reference_point( 1 ) - 0.3, 1e-13 );
}

/** @test Verify embedded exact inversion and closest-point projection. */
BOOST_AUTO_TEST_CASE( embedded_curve_closest_points )
{
    using geomap_type = GeoMap<1, 2, 2, double, Simplex>;
    auto geomap = std::make_shared<geomap_type>();
    constexpr double curvature = 0.28;
    auto mapping = []( Eigen::Matrix<double, 1, 1> const& point )
    {
        Eigen::Vector2d result;
        result << point( 0 ), curvature * point( 0 ) * point( 0 );
        return result;
    };
    auto geometry = makeGeometry( geomap, mapping );
    typename geomap_type::Inverse inverse( geomap, geometry );

    Eigen::Matrix<double, 1, 1> referencePoint;
    referencePoint << 0.23;
    auto const exactTarget = makeNode( mapping( referencePoint ) );
    auto const exact = inverse.inverse( exactTarget );
    BOOST_CHECK( exact.exact() );
    BOOST_CHECK( exact.admissible() );

    Eigen::Vector2d tangent( 1.0, 2.0 * curvature * referencePoint( 0 ) );
    Eigen::Vector2d normal( -tangent( 1 ), tangent( 0 ) );
    normal.normalize();
    Eigen::Vector2d const offCurveTarget =
        mapping( referencePoint ) + 0.12 * normal;
    auto const closest = inverse.closestPoint( makeNode( offCurveTarget ) );
    BOOST_CHECK( closest.converged() );
    BOOST_CHECK( !closest.exact() );
    BOOST_CHECK( closest.status ==
                 Feel::geometry::InverseStatus::converged_projected );
    BOOST_CHECK_SMALL(
        closest.reference_point( 0 ) - referencePoint( 0 ), 2e-10 );
    BOOST_CHECK_SMALL( closest.residual_norm - 0.12, 2e-10 );
    BOOST_CHECK_SMALL( closest.distance() - 0.12, 2e-10 );
    BOOST_CHECK( !closest.onBoundary() );
    BOOST_CHECK_SMALL( closest.gradient_norm, 2e-11 );

    Eigen::Matrix<double, 1, 1> extendedReference;
    extendedReference << 1.3;
    auto const extendedTarget = makeNode( mapping( extendedReference ) );
    auto const unrestricted = inverse.inverse( extendedTarget );
    BOOST_CHECK( unrestricted.exact() );
    BOOST_CHECK( unrestricted.location ==
                 Feel::geometry::InverseLocation::outside );

    typename geomap_type::Inverse::inverse_options_type constrained;
    constrained.warm_start = false;
    auto const boundary = inverse.closestPoint( extendedTarget, constrained );
    BOOST_CHECK( boundary.converged() );
    BOOST_CHECK( boundary.location ==
                 Feel::geometry::InverseLocation::boundary );
    BOOST_CHECK( boundary.onBoundary() );
    BOOST_CHECK_SMALL( boundary.reference_point( 0 ) - 1.0, 2e-12 );
    BOOST_CHECK_GT( boundary.residual_norm, 0.2 );

    Eigen::Matrix<double, 2, 2> projectionTargets;
    projectionTargets.col( 0 ) = offCurveTarget;
    projectionTargets.col( 1 ) = eigenMap( extendedTarget );
    auto const projectedBatch = inverse.closestPoint( projectionTargets );
    BOOST_REQUIRE_EQUAL( projectedBatch.size(), 2 );
    BOOST_CHECK( projectedBatch[0].converged() );
    BOOST_CHECK_SMALL( projectedBatch[0].distance()-0.12, 2e-10 );
    BOOST_CHECK( projectedBatch[1].converged() );
    BOOST_CHECK( projectedBatch[1].onBoundary() );
}

/** @test Verify batched nonlinear inversion and sequential warm starts. */
BOOST_AUTO_TEST_CASE( nonlinear_batch_warm_start )
{
    using geomap_type = GeoMap<2, 3, 2, double, Simplex>;
    auto geomap = std::make_shared<geomap_type>();
    auto geometry = makeGeometry(
        geomap, []( Eigen::Vector2d const& point )
        { return polynomialMap<3>( point ); } );
    typename geomap_type::Inverse inverse( geomap, geometry );

    typename geomap_type::matrix_node_t_type targets( 2, 4 );
    std::array<Eigen::Vector2d, 4> const referencePoints = {
        Eigen::Vector2d( -0.5, -0.4 ),
        Eigen::Vector2d( -0.45, -0.38 ),
        Eigen::Vector2d( -0.4, -0.36 ),
        Eigen::Vector2d( -0.35, -0.34 ) };
    for ( std::size_t point = 0; point < referencePoints.size(); ++point )
    {
        auto const physicalPoint = polynomialMap<3>( referencePoints[point] );
        targets( 0, point ) = physicalPoint( 0 );
        targets( 1, point ) = physicalPoint( 1 );
    }

    auto const results = inverse.inverse( targets );
    BOOST_REQUIRE_EQUAL( results.size(), referencePoints.size() );
    for ( std::size_t point = 0; point < results.size(); ++point )
    {
        BOOST_CHECK( results[point].exact() );
        BOOST_CHECK_SMALL(
            ( results[point].reference_point -
              referencePoints[point] )
                .norm(),
            2e-10 );
        BOOST_CHECK_EQUAL( results[point].used_warm_start, point > 0 );
    }

    auto const compatibilityCoordinates = inverse( targets );
    for ( std::size_t point = 0; point < referencePoints.size(); ++point )
        for ( int component = 0; component < 2; ++component )
            BOOST_CHECK_SMALL(
                compatibilityCoordinates( component, point ) -
                    referencePoints[point]( component ),
                2e-10 );
}

/** @test Verify explicit numerical failure diagnostics. */
BOOST_AUTO_TEST_CASE( failure_diagnostics )
{
    using geomap_type = GeoMap<2, 2, 2, double, Simplex>;
    auto geomap = std::make_shared<geomap_type>();
    auto degenerateGeometry = makeGeometry(
        geomap, []( Eigen::Vector2d const& /**/ )
        { return Eigen::Vector2d::Zero(); } );
    typename geomap_type::Inverse inverse( geomap, degenerateGeometry );

    auto const rankFailure = inverse.inverse(
        makeNode( Eigen::Vector2d( 0.4, -0.2 ) ) );
    BOOST_CHECK( rankFailure.status ==
                 Feel::geometry::InverseStatus::rank_deficient );
    BOOST_CHECK_EQUAL( rankFailure.jacobian_rank, 0 );
    BOOST_CHECK( !rankFailure.converged() );

    node_t<double> nonFiniteTarget( 2 );
    nonFiniteTarget( 0 ) = std::numeric_limits<double>::quiet_NaN();
    nonFiniteTarget( 1 ) = 0.0;
    auto const nonFinite = inverse.inverse( nonFiniteTarget );
    BOOST_CHECK( nonFinite.status ==
                 Feel::geometry::InverseStatus::non_finite );

    auto regularGeometry = makeGeometry(
        geomap, []( Eigen::Vector2d const& point )
        { return polynomialMap<2>( point ); } );
    inverse.update( regularGeometry );
    typename geomap_type::Inverse::inverse_options_type noIterations;
    noIterations.maximum_iterations = 0;
    noIterations.warm_start = false;
    auto const iterationFailure = inverse.inverse(
        makeNode( polynomialMap<2>( Eigen::Vector2d( 0.1, -0.4 ) ) ),
        noIterations );
    BOOST_CHECK( iterationFailure.status ==
                 Feel::geometry::InverseStatus::maximum_iterations );
    BOOST_CHECK( !iterationFailure.converged() );

    using evaluation_type =
        Feel::geometry::InverseMapEvaluation<double, 2, 2>;
    auto inconsistentEvaluator = []( Eigen::Vector2d const& /**/ )
    {
        return evaluation_type{ Eigen::Vector2d::Zero(),
                                Eigen::Matrix2d::Identity() };
    };
    Eigen::Vector2d const inconsistentTarget = Eigen::Vector2d::Ones();
    Eigen::Vector2d const inconsistentInitial = Eigen::Vector2d::Zero();
    auto const lineSearchFailure = Feel::geometry::solveInverseMap<2, 2>(
        inconsistentTarget, inconsistentInitial,
        inconsistentEvaluator );
    BOOST_CHECK( lineSearchFailure.status ==
                 Feel::geometry::InverseStatus::line_search_failed );

    Feel::geometry::InverseOptions<double> stagnationOptions;
    stagnationOptions.maximum_line_search_iterations = 64;
    stagnationOptions.armijo_coefficient = 1.0;
    auto const stagnationFailure = Feel::geometry::solveInverseMap<2, 2>(
        inconsistentTarget, inconsistentInitial,
        inconsistentEvaluator, stagnationOptions );
    BOOST_CHECK_MESSAGE(
        stagnationFailure.status == Feel::geometry::InverseStatus::stagnated,
        "expected stagnation, got " << Feel::geometry::inverseStatusName( stagnationFailure.status ) );
}

BOOST_AUTO_TEST_SUITE_END()
