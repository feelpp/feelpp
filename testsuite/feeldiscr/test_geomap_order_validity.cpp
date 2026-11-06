/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 * @file test_geomap_order_validity.cpp
 * @brief Analytic runtime-order, mapping-quality, and validity tests for GeoMap.
 */

#define BOOST_TEST_MODULE geomap order and validity testsuite
#include <feel/feelcore/testsuite.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmshenums.hpp>
#include <feel/feelmesh/geond.hpp>
#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/geometrydescriptor.hpp>
#include <feel/feelpoly/geometrypointtopology.hpp>
#include <feel/feelpoly/geomap.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

using namespace Feel;

namespace
{
/**
 * @brief Evaluate a known polynomial physical mapping used as a test oracle.
 * @tparam Dim reference dimension
 * @tparam RealDim ambient physical dimension
 * @tparam ReferencePoint reference-point vector type
 * @param point reference point
 * @param order polynomial order to retain in the oracle
 * @return physical point produced by the analytic mapping
 */
template <int Dim, int RealDim, typename ReferencePoint>
auto analyticMap( ReferencePoint const& point, int order )
{
    Eigen::Matrix<double, RealDim, 1> result;
    const bool quadratic = order >= 2;
    const bool cubic = order >= 3;

    if constexpr ( Dim == 1 && RealDim == 1 )
    {
        const double x = point( 0 );
        result( 0 ) = 0.3 + 1.2 * x +
                      ( quadratic ? 0.08 * x * x : 0.0 ) +
                      ( cubic ? 0.01 * x * x * x : 0.0 );
    }
    else if constexpr ( Dim == 1 && RealDim == 2 )
    {
        const double x = point( 0 );
        result( 0 ) = 0.3 + 1.2 * x + ( cubic ? 0.01 * x * x * x : 0.0 );
        result( 1 ) = -0.2 + 0.4 * x + ( quadratic ? 0.1 * x * x : 0.0 );
    }
    else if constexpr ( Dim == 1 && RealDim == 3 )
    {
        const double x = point( 0 );
        result( 0 ) = 0.1 + x + ( cubic ? 0.015*x*x*x : 0.0 );
        result( 1 ) = -0.2 + 0.5*x + ( quadratic ? 0.08*x*x : 0.0 );
        result( 2 ) = 0.3 - 0.25*x + ( quadratic ? 0.04*x*x : 0.0 ) +
                      ( cubic ? 0.01*x*x*x : 0.0 );
    }
    else if constexpr ( Dim == 2 && RealDim == 2 )
    {
        const double x = point( 0 );
        const double y = point( 1 );
        result( 0 ) = 0.2 + 1.1 * x + 0.1 * y +
                      ( quadratic ? 0.04 * x * y + 0.03 * y * y : 0.0 ) +
                      ( cubic ? 0.006 * x * x * x : 0.0 );
        result( 1 ) = -0.1 + 0.05 * x + 0.9 * y +
                      ( quadratic ? 0.025 * x * x - 0.02 * x * y : 0.0 ) +
                      ( cubic ? 0.005 * y * y * y : 0.0 );
    }
    else if constexpr ( Dim == 2 && RealDim == 3 )
    {
        const double x = point( 0 );
        const double y = point( 1 );
        result( 0 ) = x + ( quadratic ? 0.03 * x * y : 0.0 ) +
                      ( cubic ? 0.004 * x * x * x : 0.0 );
        result( 1 ) = y + ( quadratic ? 0.02 * x * x : 0.0 ) +
                      ( cubic ? 0.004 * y * y * y : 0.0 );
        result( 2 ) = 0.2 + 0.1 * x - 0.05 * y +
                      ( quadratic ? 0.08 * ( x * x + y * y ) : 0.0 );
    }
    else if constexpr ( Dim == 3 && RealDim == 3 )
    {
        const double x = point( 0 );
        const double y = point( 1 );
        const double z = point( 2 );
        result( 0 ) = 0.1 + 1.1 * x + 0.05 * y + 0.02 * z +
                      ( quadratic ? 0.02 * x * y : 0.0 ) +
                      ( cubic ? 0.003 * x * x * x : 0.0 );
        result( 1 ) = -0.2 + 0.03 * x + 0.95 * y + 0.04 * z +
                      ( quadratic ? 0.015 * y * z : 0.0 ) +
                      ( cubic ? 0.003 * y * y * y : 0.0 );
        result( 2 ) = 0.3 + 0.02 * x - 0.01 * y + 1.05 * z +
                      ( quadratic ? 0.01 * x * z : 0.0 ) +
                      ( cubic ? 0.003 * z * z * z : 0.0 );
    }
    return result;
}

/**
 * @brief Evaluate the exact Jacobian of `analyticMap`.
 * @tparam Dim reference dimension
 * @tparam RealDim ambient physical dimension
 * @tparam ReferencePoint reference-point vector type
 * @param point reference point
 * @param order polynomial order to retain in the oracle
 * @return exact `RealDim`-by-`Dim` physical Jacobian matrix
 */
template <int Dim, int RealDim, typename ReferencePoint>
auto analyticJacobian( ReferencePoint const& point, int order )
{
    Eigen::Matrix<double, RealDim, Dim> result;
    result.setZero();
    const bool quadratic = order >= 2;
    const bool cubic = order >= 3;

    if constexpr ( Dim == 1 && RealDim == 1 )
    {
        const double x = point( 0 );
        result( 0, 0 ) = 1.2 + ( quadratic ? 0.16 * x : 0.0 ) +
                         ( cubic ? 0.03 * x * x : 0.0 );
    }
    else if constexpr ( Dim == 1 && RealDim == 2 )
    {
        const double x = point( 0 );
        result( 0, 0 ) = 1.2 + ( cubic ? 0.03 * x * x : 0.0 );
        result( 1, 0 ) = 0.4 + ( quadratic ? 0.2 * x : 0.0 );
    }
    else if constexpr ( Dim == 1 && RealDim == 3 )
    {
        const double x = point( 0 );
        result( 0, 0 ) = 1.0 + ( cubic ? 0.045*x*x : 0.0 );
        result( 1, 0 ) = 0.5 + ( quadratic ? 0.16*x : 0.0 );
        result( 2, 0 ) = -0.25 + ( quadratic ? 0.08*x : 0.0 ) +
                         ( cubic ? 0.03*x*x : 0.0 );
    }
    else if constexpr ( Dim == 2 && RealDim == 2 )
    {
        const double x = point( 0 );
        const double y = point( 1 );
        result( 0, 0 ) = 1.1 + ( quadratic ? 0.04 * y : 0.0 ) +
                         ( cubic ? 0.018 * x * x : 0.0 );
        result( 0, 1 ) = 0.1 + ( quadratic ? 0.04 * x + 0.06 * y : 0.0 );
        result( 1, 0 ) = 0.05 + ( quadratic ? 0.05 * x - 0.02 * y : 0.0 );
        result( 1, 1 ) = 0.9 - ( quadratic ? 0.02 * x : 0.0 ) +
                         ( cubic ? 0.015 * y * y : 0.0 );
    }
    else if constexpr ( Dim == 2 && RealDim == 3 )
    {
        const double x = point( 0 );
        const double y = point( 1 );
        result( 0, 0 ) = 1.0 + ( quadratic ? 0.03 * y : 0.0 ) +
                         ( cubic ? 0.012 * x * x : 0.0 );
        result( 0, 1 ) = quadratic ? 0.03 * x : 0.0;
        result( 1, 0 ) = quadratic ? 0.04 * x : 0.0;
        result( 1, 1 ) = 1.0 + ( cubic ? 0.012 * y * y : 0.0 );
        result( 2, 0 ) = 0.1 + ( quadratic ? 0.16 * x : 0.0 );
        result( 2, 1 ) = -0.05 + ( quadratic ? 0.16 * y : 0.0 );
    }
    else if constexpr ( Dim == 3 && RealDim == 3 )
    {
        const double x = point( 0 );
        const double y = point( 1 );
        const double z = point( 2 );
        result << 1.1 + ( quadratic ? 0.02 * y : 0.0 ) + ( cubic ? 0.009 * x * x : 0.0 ),
            0.05 + ( quadratic ? 0.02 * x : 0.0 ),
            0.02,
            0.03,
            0.95 + ( quadratic ? 0.015 * z : 0.0 ) + ( cubic ? 0.009 * y * y : 0.0 ),
            0.04 + ( quadratic ? 0.015 * y : 0.0 ),
            0.02 + ( quadratic ? 0.01 * z : 0.0 ),
            -0.01,
            1.05 + ( quadratic ? 0.01 * x : 0.0 ) + ( cubic ? 0.009 * z * z : 0.0 );
    }
    return result;
}

/**
 * @brief Check analytic correctness and static/runtime parity for one map.
 *
 * The check covers mapped points, the physical Jacobian, integration and
 * signed Jacobians, inverse-transpose identity, singular values, condition
 * number, validity, and propagation of runtime order to the inverse map.
 *
 * @tparam Dim reference dimension
 * @tparam RealDim ambient physical dimension
 * @tparam Order polynomial mapping order
 * @tparam Convex reference-entity family
 */
template <int Dim, int RealDim, int Order,
          template <int, int, int> class Convex>
void checkAnalyticStaticRuntimeParity()
{
    using static_geomap_type = GeoMap<Dim, Order, RealDim, double, Convex>;
    using runtime_geomap_type = GeoMap<Dim, Dynamic, RealDim, double, Convex>;
    using shape_type = Convex<Dim, Order, RealDim>;
    using element_type = GeoND<Dim, shape_type>;
    using point_type = typename element_type::point_type;
    using node_type = typename point_type::node_type;

    constexpr size_type contextFlags =
        vm::POINT | vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE;
    constexpr double tolerance = 5e-11;

    auto staticGeomap = std::make_shared<static_geomap_type>();
    auto runtimeGeomap = std::make_shared<runtime_geomap_type>( RuntimeOrder{ Order } );
    BOOST_REQUIRE_EQUAL( runtimeGeomap->order(), Order );
    BOOST_REQUIRE_EQUAL( runtimeGeomap->nbPoints(), staticGeomap->nbPoints() );

    element_type element;
    std::vector<point_type> pointStorage;
    pointStorage.reserve( staticGeomap->nbPoints() );
    for ( uint16_type i = 0; i < staticGeomap->nbPoints(); ++i )
    {
        auto const referenceNode = staticGeomap->refNode( i );
        auto const physicalNode = analyticMap<Dim, RealDim>( referenceNode, Order );
        node_type node( RealDim );
        for ( int component = 0; component < RealDim; ++component )
            node( component ) = physicalNode( component );
        pointStorage.emplace_back( static_cast<uint32_type>( i ), node, false,
                                   i < shape_type::numVertices );
        element.setPoint( i, pointStorage.back() );
    }

    auto const referencePoints = staticGeomap->points();
    auto staticPrecompute = staticGeomap->preCompute( staticGeomap, referencePoints );
    auto runtimePrecompute = runtimeGeomap->preCompute( runtimeGeomap, referencePoints );
    auto staticContext = staticGeomap->template context<contextFlags>( element, staticPrecompute );
    auto runtimeContext = runtimeGeomap->template context<contextFlags>( element, runtimePrecompute );

    static_assert( GeometricMappingContextConcept<decltype( *staticContext )> );
    static_assert( GeometricMappingContextConcept<decltype( *runtimeContext )> );
    if constexpr ( Dim < RealDim )
    {
        static_assert( EmbeddedGeometricMappingContextConcept<
                       decltype( *staticContext )> );
        static_assert( EmbeddedGeometricMappingContextConcept<
                       decltype( *runtimeContext )> );
    }
    BOOST_REQUIRE( staticContext->isValid() );
    BOOST_REQUIRE( runtimeContext->isValid() );

    if constexpr ( static_geomap_type::is_linear )
    {
        BOOST_CHECK( staticContext->usesAffineStorage() );
        BOOST_CHECK( runtimeContext->usesAffineStorage() );
        BOOST_CHECK_EQUAL( staticContext->nComputedPoints(), 1 );
        BOOST_CHECK_EQUAL( runtimeContext->nComputedPoints(), 1 );
        BOOST_CHECK_EQUAL( runtimeContext->jacobianMatrices().size(), 1 );
        BOOST_CHECK_EQUAL(
            runtimeContext->inverseJacobianTransposeMatrices().size(), 1 );
        BOOST_CHECK_EQUAL(
            runtimeContext->integrationJacobians().size(), 1 );
        for ( int q = 0; q < runtimeContext->nPoints(); ++q )
        {
            BOOST_CHECK_EQUAL(
                runtimeContext->jacobianMatrix( q ).data(),
                runtimeContext->jacobianMatrix( 0 ).data() );
            for ( int component = 0; component < RealDim; ++component )
                for ( int first = 0; first < Dim; ++first )
                    for ( int second = 0; second < Dim; ++second )
                        BOOST_CHECK_SMALL(
                            runtimeContext->hessian( q )(
                                component, first, second ),
                            tolerance );
        }
    }
    else
    {
        BOOST_CHECK( !runtimeContext->usesAffineStorage() );
        BOOST_CHECK_EQUAL( runtimeContext->nComputedPoints(),
                           runtimeContext->nPoints() );

        runtimeContext->setGeometryEvaluationStrategy(
            geometry::geometry_evaluation_strategy::pointwise );
        runtimeContext->template update<contextFlags>( element );
        std::vector<typename runtime_geomap_type::template Context<
            element_type>::eigen_matrix_np_type>
            pointwiseK;
        std::vector<typename runtime_geomap_type::template Context<
            element_type>::eigen_matrix_np_type>
            pointwiseB;
        std::vector<double> pointwiseJ;
        std::vector<typename runtime_geomap_type::template Context<
            element_type>::hessian_type>
            pointwiseHessian;
        for ( int q = 0; q < runtimeContext->nPoints(); ++q )
        {
            pointwiseK.push_back( runtimeContext->jacobianMatrix( q ) );
            pointwiseB.push_back(
                runtimeContext->inverseJacobianTranspose( q ) );
            pointwiseJ.push_back(
                runtimeContext->integrationJacobian( q ) );
            pointwiseHessian.push_back( runtimeContext->hessian( q ) );
        }
        auto const* jacobianStorage =
            runtimeContext->jacobianMatrices().data();
        auto const* inverseStorage =
            runtimeContext->inverseJacobianTransposeMatrices().data();
        auto const* measureStorage =
            runtimeContext->integrationJacobians().data();

        runtimeContext->setGeometryEvaluationStrategy(
            geometry::geometry_evaluation_strategy::batched );
        runtimeContext->template update<contextFlags>( element );
        BOOST_CHECK_EQUAL( runtimeContext->jacobianMatrices().data(),
                           jacobianStorage );
        BOOST_CHECK_EQUAL(
            runtimeContext->inverseJacobianTransposeMatrices().data(),
            inverseStorage );
        BOOST_CHECK_EQUAL( runtimeContext->integrationJacobians().data(),
                           measureStorage );
        for ( int q = 0; q < runtimeContext->nPoints(); ++q )
        {
            BOOST_CHECK_SMALL(
                ( runtimeContext->jacobianMatrix( q ) - pointwiseK[q] )
                    .norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                ( runtimeContext->inverseJacobianTranspose( q ) -
                  pointwiseB[q] )
                    .norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                runtimeContext->integrationJacobian( q ) - pointwiseJ[q],
                tolerance );
            for ( int component = 0; component < RealDim; ++component )
                for ( int first = 0; first < Dim; ++first )
                    for ( int second = 0; second < Dim; ++second )
                        BOOST_CHECK_SMALL(
                            runtimeContext->hessian( q )(
                                component, first, second ) -
                                pointwiseHessian[q](
                                    component, first, second ),
                            tolerance );
        }
    }

    double expectedMinimumJacobian = std::numeric_limits<double>::infinity();
    for ( int q = 0; q < staticContext->nPoints(); ++q )
    {
        auto const referencePoint = ublas::column( referencePoints, q );
        auto const expectedPoint = analyticMap<Dim, RealDim>( referencePoint, Order );
        auto const expectedJacobian = analyticJacobian<Dim, RealDim>( referencePoint, Order );
        const double expectedMeasure = std::sqrt(
            ( expectedJacobian.transpose() * expectedJacobian ).determinant() );
        expectedMinimumJacobian = std::min( expectedMinimumJacobian, expectedMeasure );

        for ( int component = 0; component < RealDim; ++component )
        {
            BOOST_CHECK_SMALL( staticContext->xReal( q )( component ) -
                                   expectedPoint( component ),
                               tolerance );
            BOOST_CHECK_SMALL( staticContext->xReal( q )( component ) -
                                   runtimeContext->xReal( q )( component ),
                               tolerance );
        }
        BOOST_CHECK_SMALL( ( staticContext->jacobianMatrix( q ) - expectedJacobian ).norm(),
                           tolerance );
        BOOST_CHECK_SMALL( ( staticContext->jacobianMatrix( q ) -
                             runtimeContext->jacobianMatrix( q ) )
                               .norm(),
                           tolerance );
        BOOST_CHECK_SMALL( staticContext->integrationJacobian( q ) - expectedMeasure,
                           tolerance );
        BOOST_CHECK_SMALL( staticContext->integrationJacobian( q ) -
                               runtimeContext->integrationJacobian( q ),
                           tolerance );
        BOOST_CHECK_SMALL(
            ( staticContext->jacobianMatrix( q ).transpose() *
                  staticContext->inverseJacobianTranspose( q ) -
              Eigen::Matrix<double, Dim, Dim>::Identity() )
                .norm(),
            tolerance );

        if constexpr ( Dim < RealDim )
        {
            auto const staticMetric = staticContext->metricTensor( q );
            auto const staticInverseMetric =
                staticContext->inverseMetricTensor( q );
            auto const staticTangent = staticContext->tangentProjector( q );
            auto const staticNormal = staticContext->normalProjector( q );
            BOOST_CHECK_SMALL(
                ( staticMetric-staticContext->jacobianMatrix( q ).transpose()*
                                    staticContext->jacobianMatrix( q ) ).norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                ( staticMetric*staticInverseMetric-
                  Eigen::Matrix<double, Dim, Dim>::Identity() ).norm(),
                tolerance );
            BOOST_CHECK_SMALL( ( staticTangent*staticTangent-staticTangent ).norm(),
                               tolerance );
            BOOST_CHECK_SMALL( ( staticNormal*staticNormal-staticNormal ).norm(),
                               tolerance );
            BOOST_CHECK_SMALL( ( staticTangent*staticNormal ).norm(), tolerance );
            BOOST_CHECK_SMALL(
                ( staticTangent+staticNormal-
                  Eigen::Matrix<double, RealDim, RealDim>::Identity() ).norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                staticContext->manifoldJacobian( q )-
                    staticContext->integrationJacobian( q ),
                tolerance );
            BOOST_CHECK_SMALL(
                ( staticMetric-runtimeContext->metricTensor( q ) ).norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                ( staticInverseMetric-
                  runtimeContext->inverseMetricTensor( q ) ).norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                ( staticTangent-runtimeContext->tangentProjector( q ) ).norm(),
                tolerance );
            BOOST_CHECK_SMALL(
                ( staticNormal-runtimeContext->normalProjector( q ) ).norm(),
                tolerance );

            auto const staticNormalFrame = staticContext->normalFrame( q );
            BOOST_CHECK_SMALL(
                ( staticContext->jacobianMatrix( q ).transpose()*
                  staticNormalFrame ).norm(), tolerance );
            BOOST_CHECK_SMALL(
                ( staticNormalFrame*staticNormalFrame.transpose()-
                  staticNormal ).norm(), tolerance );

            if constexpr ( RealDim == Dim+1 && RealDim <= 3 )
            {
                BOOST_CHECK_SMALL(
                    ( staticContext->unitNormal( q )-
                      runtimeContext->unitNormal( q ) ).norm(), tolerance );
                BOOST_CHECK_SMALL(
                    ( staticContext->secondFundamentalForm( q )-
                      runtimeContext->secondFundamentalForm( q ) ).norm(),
                    tolerance );
                BOOST_CHECK_SMALL(
                    ( staticContext->shapeOperator( q )-
                      runtimeContext->shapeOperator( q ) ).norm(), tolerance );
                BOOST_CHECK_SMALL(
                    staticContext->meanCurvature( q )-
                        runtimeContext->meanCurvature( q ), tolerance );
                if constexpr ( Dim == 2 && RealDim == 3 )
                    BOOST_CHECK_SMALL(
                        staticContext->gaussianCurvature( q )-
                            runtimeContext->gaussianCurvature( q ), tolerance );
            }
        }

        auto const staticQuality = staticContext->mappingQuality( q );
        auto const runtimeQuality = runtimeContext->mappingQuality( q );
        BOOST_CHECK( staticQuality.valid );
        BOOST_CHECK( staticQuality.finite );
        BOOST_CHECK( staticQuality.orientation_preserving );
        BOOST_CHECK_CLOSE_FRACTION( staticQuality.condition_number,
                                    runtimeQuality.condition_number, tolerance );
        BOOST_CHECK_CLOSE_FRACTION( staticQuality.minimum_singular_value,
                                    runtimeQuality.minimum_singular_value, tolerance );
        if constexpr ( Dim == RealDim )
            BOOST_CHECK_CLOSE_FRACTION( staticContext->signedJacobian( q ),
                                        expectedJacobian.determinant(), tolerance );
        else
            BOOST_CHECK_CLOSE_FRACTION( staticContext->signedJacobian( q ),
                                        expectedMeasure, tolerance );
    }
    BOOST_CHECK_CLOSE_FRACTION( staticContext->minimumJacobian(),
                                expectedMinimumJacobian, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( staticContext->minimumJacobian(),
                                runtimeContext->minimumJacobian(), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( staticContext->maximumConditionNumber(),
                                runtimeContext->maximumConditionNumber(), tolerance );

    typename runtime_geomap_type::Inverse runtimeInverse( runtimeGeomap, element );
    BOOST_CHECK_EQUAL( runtimeInverse.geometricMapping()->order(), Order );
}

/**
 * @brief Validate a runtime-only map above the static build envelope.
 *
 * This helper instantiates only `GeoMap<..., Dynamic>` and a dynamic-order
 * entity. It checks that runtime tabulation, Context derivatives, quality,
 * and Inverse construction remain mathematical at an order for which the
 * build provides no corresponding static mesh instantiation.
 *
 * @tparam Dim reference dimension
 * @tparam RealDim ambient physical dimension
 * @tparam Convex reference-entity family
 * @param order runtime geometric polynomial order
 */
template <int Dim, int RealDim,
          template <int, int, int> class Convex>
void checkAnalyticRuntimeOrder( uint16_type order )
{
    using geomap_type =
        GeoMap<Dim, Dynamic, RealDim, double, Convex>;
    using shape_type = Convex<Dim, Dynamic, RealDim>;
    using element_type = GeoND<Dim, shape_type>;
    using point_type = typename element_type::point_type;
    using node_type = typename point_type::node_type;

    constexpr size_type contextFlags =
        vm::POINT | vm::JACOBIAN | vm::KB |
        vm::SECOND_DERIVATIVE;
    constexpr double tolerance = 2e-9;

    auto geomap =
        std::make_shared<geomap_type>(
            RuntimeOrder{ order } );
    element_type element{ RuntimeOrder{ order } };
    std::vector<point_type> pointStorage;
    pointStorage.reserve( geomap->nbPoints() );
    for ( uint16_type point = 0;
          point < geomap->nbPoints(); ++point )
    {
        auto const physicalNode =
            analyticMap<Dim, RealDim>(
                geomap->refNode( point ), order );
        node_type node( RealDim );
        for ( int component = 0;
              component < RealDim; ++component )
            node( component ) = physicalNode( component );
        pointStorage.emplace_back(
            static_cast<uint32_type>( point ), node, false,
            point < shape_type::numVertices );
        element.setPoint( point, pointStorage.back() );
    }

    auto const referencePoints = geomap->points();
    auto precompute =
        geomap->preCompute( geomap, referencePoints );
    auto context =
        geomap->template context<contextFlags>(
            element, precompute );

    BOOST_REQUIRE_EQUAL( geomap->order(), order );
    BOOST_REQUIRE( context->isValid() );
    BOOST_CHECK( !context->usesAffineStorage() );
    BOOST_CHECK_EQUAL(
        context->nComputedPoints(), context->nPoints() );

    for ( int q = 0; q < context->nPoints(); ++q )
    {
        auto const referencePoint =
            ublas::column( referencePoints, q );
        auto const expectedPoint =
            analyticMap<Dim, RealDim>(
                referencePoint, order );
        auto const expectedJacobian =
            analyticJacobian<Dim, RealDim>(
                referencePoint, order );
        for ( int component = 0;
              component < RealDim; ++component )
            BOOST_CHECK_SMALL(
                context->xReal( q )( component ) -
                    expectedPoint( component ),
                tolerance );
        BOOST_CHECK_SMALL(
            ( context->jacobianMatrix( q ) -
              expectedJacobian )
                .norm(),
            tolerance );
        BOOST_CHECK( context->mappingQuality( q ).valid );
    }

    typename geomap_type::Inverse inverse{
        geomap, element };
    BOOST_CHECK_EQUAL(
        inverse.geometricMapping()->order(), order );
}

/**
 * @brief Build a two-dimensional simplex context from a supplied mapping.
 * @tparam Order compile-time geometric order
 * @tparam Mapping callable mapping reference points to physical points
 * @tparam Check callable receiving the fully updated geometric context
 * @param mapping physical mapping used to place all geometric nodes
 * @param check assertions executed on the resulting context
 */
template <int Order, typename Mapping, typename Check>
void withSimplex2DContext( Mapping&& mapping, Check&& check )
{
    using geomap_type = GeoMap<2, Order, 2>;
    using shape_type = Simplex<2, Order, 2>;
    using element_type = GeoND<2, shape_type>;
    using point_type = typename element_type::point_type;
    using node_type = typename point_type::node_type;
    constexpr size_type contextFlags = vm::POINT | vm::JACOBIAN | vm::KB;

    auto geomap = std::make_shared<geomap_type>();
    element_type element;
    std::vector<point_type> pointStorage;
    pointStorage.reserve( geomap->nbPoints() );
    for ( uint16_type i = 0; i < geomap->nbPoints(); ++i )
    {
        auto const physicalPoint = mapping( geomap->refNode( i ) );
        node_type node( 2 );
        node( 0 ) = physicalPoint( 0 );
        node( 1 ) = physicalPoint( 1 );
        pointStorage.emplace_back( static_cast<uint32_type>( i ), node, false,
                                   i < shape_type::numVertices );
        element.setPoint( i, pointStorage.back() );
    }
    auto const points = geomap->points();
    auto precompute = geomap->preCompute( geomap, points );
    auto context = geomap->template context<contextFlags>( element, precompute );
    check( *context );
}

/**
 * @brief Evaluate one embedded triangular surface patch with static and runtime order.
 * @tparam Order polynomial geometric order
 * @tparam Mapping analytic surface parametrization sampled at geometric nodes
 * @tparam Check callable receiving static and runtime contexts at the barycenter
 * @param mapping physical surface parametrization
 * @param check assertions or error extraction callback
 */
template <int Order, typename Mapping, typename Check>
void withEmbeddedSurfaceContext( Mapping&& mapping, Check&& check )
{
    using static_geomap_type = GeoMap<2, Order, 3>;
    using runtime_geomap_type = GeoMap<2, Dynamic, 3>;
    using shape_type = Simplex<2, Order, 3>;
    using element_type = GeoND<2, shape_type>;
    using point_type = typename element_type::point_type;
    using node_type = typename point_type::node_type;
    constexpr size_type contextFlags =
        vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE;

    auto staticGeomap = std::make_shared<static_geomap_type>();
    auto runtimeGeomap =
        std::make_shared<runtime_geomap_type>( RuntimeOrder{ Order } );
    element_type element;
    std::vector<point_type> pointStorage;
    pointStorage.reserve( staticGeomap->nbPoints() );
    for ( uint16_type i = 0; i < staticGeomap->nbPoints(); ++i )
    {
        auto const physicalPoint = mapping( staticGeomap->refNode( i ) );
        node_type node( 3 );
        for ( int component = 0; component < 3; ++component )
            node( component ) = physicalPoint( component );
        pointStorage.emplace_back( static_cast<uint32_type>( i ), node, false,
                                   i < shape_type::numVertices );
        element.setPoint( i, pointStorage.back() );
    }

    typename static_geomap_type::matrix_node_t_type evaluationPoint( 2, 1 );
    evaluationPoint( 0, 0 ) = -1.0/3.0;
    evaluationPoint( 1, 0 ) = -1.0/3.0;
    auto staticPrecompute = staticGeomap->preCompute(
        staticGeomap, evaluationPoint );
    auto runtimePrecompute = runtimeGeomap->preCompute(
        runtimeGeomap, evaluationPoint );
    auto staticContext = staticGeomap->template context<contextFlags>(
        element, staticPrecompute );
    auto runtimeContext = runtimeGeomap->template context<contextFlags>(
        element, runtimePrecompute );
    check( *staticContext, *runtimeContext );
}

/** @brief Curvature errors of one polynomial surface patch. */
struct curvature_error_type
{
    double mean = 0.0;     /**< Absolute signed-mean-curvature error. */
    double gaussian = 0.0; /**< Absolute Gaussian-curvature error. */
};

/**
 * @brief Evaluate static/runtime curvature and return the static approximation error.
 * @tparam Order polynomial geometric order
 * @param mapping analytic parametrization used for geometric interpolation
 * @param exactMean exact signed mean curvature at the evaluation point
 * @param exactGaussian exact Gaussian curvature at the evaluation point
 */
template <int Order, typename Mapping>
curvature_error_type surfaceCurvatureError(
    Mapping&& mapping, double exactMean, double exactGaussian )
{
    curvature_error_type error;
    withEmbeddedSurfaceContext<Order>(
        std::forward<Mapping>( mapping ),
        [&]( auto const& staticContext, auto const& runtimeContext )
        {
            constexpr double parityTolerance = 2e-10;
            BOOST_CHECK_SMALL( staticContext.meanCurvature( 0 )-
                                   runtimeContext.meanCurvature( 0 ),
                               parityTolerance );
            BOOST_CHECK_SMALL( staticContext.gaussianCurvature( 0 )-
                                   runtimeContext.gaussianCurvature( 0 ),
                               parityTolerance );
            error.mean = std::abs( staticContext.meanCurvature( 0 )-
                                   exactMean );
            error.gaussian = std::abs(
                staticContext.gaussianCurvature( 0 )-exactGaussian );
        } );
    return error;
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( geomap_order_validity_suite )

/**
 * @test Verify the affine default and invalid runtime-order rejection.
 */
BOOST_AUTO_TEST_CASE( default_and_runtime_order_contract )
{
    using default_geomap_type = GeoMap<2>;
    static_assert( std::same_as<default_geomap_type, GeoMap<2, 1, 2>> );
    static_assert( default_geomap_type::is_order_static );
    static_assert( default_geomap_type::is_linear );
    default_geomap_type defaultGeomap;
    BOOST_CHECK_EQUAL( defaultGeomap.order(), 1 );

    BOOST_CHECK_THROW( RuntimeOrder{ -1 }, std::invalid_argument );
    BOOST_CHECK_THROW( ( GeoMap<2, Dynamic>{ RuntimeOrder{ 0 } } ), std::invalid_argument );

    using dynamic_shape_type = Simplex<2, Dynamic>;
    using dynamic_entity_type = GeoND<2, dynamic_shape_type>;
    using static_entity_type = GeoND<2, Simplex<2, 2>>;
    static_assert( dynamic_entity_type::is_order_dynamic );
    static_assert( dynamic_entity_type::gm_type::is_order_dynamic );
    static_assert( HasPointCardinality<dynamic_entity_type> );
    static_assert( HasPointCardinality<static_entity_type> );

    dynamic_entity_type dynamicEntity( RuntimeOrder{ 6 } );
    dynamic_entity_type parityDynamicEntity( RuntimeOrder{ 2 } );
    static_entity_type staticEntity;
    BOOST_CHECK_EQUAL( dynamicEntity.order(), 6 );
    BOOST_CHECK_EQUAL( dynamicEntity.nPoints(), 28 );
    BOOST_CHECK_EQUAL(
        staticEntity.order(), parityDynamicEntity.order() );
    BOOST_CHECK_EQUAL(
        staticEntity.nPoints(), parityDynamicEntity.nPoints() );

    GeoEntity<Simplex<2, Dynamic>> dynamicBaseEntity(
        0, RuntimeOrder{ 4 } );
    BOOST_CHECK_EQUAL( dynamicBaseEntity.order(), 4 );
    BOOST_CHECK_EQUAL( dynamicBaseEntity.nPoints(), 15 );

    GeoElement0D<2> pointEntity;
    BOOST_CHECK_EQUAL( pointEntity.nPoints(), 1 );

    using runtime_descriptor_type =
        geometry::GeometryDescriptor<
            2, Dynamic,
            geometry::GeometryCellKind::Simplex>;
    static_assert( std::same_as<
                   runtime_descriptor_type::local_coefficient_index_type,
                   uint16_type> );
    static_assert(
        runtime_descriptor_type::max_local_coefficient_count ==
        Feel::detail::maxLocalIndexCardinality );
    runtime_descriptor_type runtimeDescriptor{
        RuntimeOrder{ 6 } };
    BOOST_CHECK_EQUAL( runtimeDescriptor.order(), 6 );
    BOOST_CHECK_EQUAL(
        runtimeDescriptor.coefficientCount(), 28 );
    BOOST_CHECK_EQUAL(
        runtimeDescriptor.measureQuadratureOrder(), 10 );
    BOOST_CHECK_THROW(
        runtime_descriptor_type( RuntimeOrder{ 0 } ),
        std::invalid_argument );

    using runtime_mesh_type =
        Mesh<Simplex<2, Dynamic>>;
    auto runtimeMesh =
        runtime_mesh_type::New(
            RuntimeOrder{ 6 },
            geometry::GeometryRepresentationKind::BernsteinBezier );
    static_assert( runtime_mesh_type::is_order_dynamic );
    BOOST_REQUIRE( runtimeMesh );
    BOOST_CHECK_EQUAL( runtimeMesh->order(), 6 );
    BOOST_CHECK(
        runtimeMesh->geometryRepresentation() ==
        geometry::GeometryRepresentationKind::BernsteinBezier );
    BOOST_CHECK_EQUAL(
        runtimeMesh->geometryDescriptor().coefficientCount(),
        28 );
    BOOST_REQUIRE( runtimeMesh->gm() );
    BOOST_CHECK_EQUAL( runtimeMesh->gm()->order(), 6 );
    BOOST_CHECK_EQUAL(
        runtimeMesh->geometryPointTopology().pointsOnFace(),
        7 );

    GeoNDCommon<dynamic_entity_type> commonData{
        runtimeMesh.get() };
    dynamic_entity_type attachedEntity{
        RuntimeOrder{ 3 } };
    BOOST_CHECK_EQUAL( attachedEntity.order(), 3 );
    BOOST_CHECK_EQUAL( attachedEntity.nPoints(), 10 );
    attachedEntity.setCommonData( &commonData );
    BOOST_CHECK_EQUAL( attachedEntity.order(), runtimeMesh->order() );
    BOOST_CHECK_EQUAL( attachedEntity.nPoints(), 28 );
    for ( uint16_type face = 0;
          face < runtime_mesh_type::element_type::numTopologicalFaces;
          ++face )
    {
        std::vector<bool> seen( 28, false );
        for ( uint16_type point = 0; point < 7; ++point )
        {
            auto const cellPoint =
                runtimeMesh->geometryPointTopology()
                    .faceToPoint( face, point );
            BOOST_REQUIRE_LT( cellPoint, 28 );
            BOOST_CHECK( !seen[cellPoint] );
            seen[cellPoint] = true;
        }
    }

    auto dynamicElement = runtimeMesh->makeElement();
    auto dynamicFace = runtimeMesh->makeFace();
    BOOST_CHECK_EQUAL( dynamicElement.order(), 6 );
    BOOST_CHECK_EQUAL( dynamicElement.nPoints(), 28 );
    BOOST_CHECK_EQUAL( dynamicFace.order(), 6 );
    BOOST_CHECK_EQUAL( dynamicFace.nPoints(), 7 );
    Feel::detail::ReversePoint<
        Simplex<2, Dynamic>> reversePoint;
    BOOST_CHECK_EQUAL(
        reversePoint.operate( 0, 28 ), 2 );
    BOOST_CHECK_EQUAL(
        reversePoint.operate( 2, 28 ), 0 );
    BOOST_CHECK_EQUAL(
        reversePoint.operate( 3, 28 ), 27 );
    BOOST_CHECK_EQUAL(
        reversePoint.operate( 27, 28 ), 3 );

    using runtime_mesh_3d_type =
        Mesh<Simplex<3, Dynamic>>;
    auto runtimeMesh3d =
        runtime_mesh_3d_type::New( RuntimeOrder{ 6 } );
    auto dynamicElement3d = runtimeMesh3d->makeElement();
    auto dynamicFace3d = runtimeMesh3d->makeFace();
    auto dynamicEdge3d = runtimeMesh3d->makeEdge();
    BOOST_CHECK_EQUAL( dynamicElement3d.order(), 6 );
    BOOST_CHECK_EQUAL( dynamicElement3d.nPoints(), 84 );
    BOOST_CHECK_EQUAL( dynamicFace3d.order(), 6 );
    BOOST_CHECK_EQUAL( dynamicFace3d.nPoints(), 28 );
    BOOST_CHECK_EQUAL( dynamicEdge3d.order(), 6 );
    BOOST_CHECK_EQUAL( dynamicEdge3d.nPoints(), 7 );
    BOOST_CHECK_EQUAL(
        runtimeMesh3d->geometryPointTopology().pointsOnFace(),
        28 );
    BOOST_CHECK_EQUAL(
        runtimeMesh3d->geometryPointTopology().pointsOnEdge(),
        7 );
    for ( uint16_type face = 0;
          face < runtime_mesh_3d_type::element_type::numTopologicalFaces;
          ++face )
    {
        std::vector<bool> seen( 84, false );
        for ( uint16_type point = 0; point < 28; ++point )
        {
            auto const cellPoint =
                runtimeMesh3d->geometryPointTopology()
                    .faceToPoint( face, point );
            BOOST_REQUIRE_LT( cellPoint, 84 );
            BOOST_CHECK( !seen[cellPoint] );
            seen[cellPoint] = true;
        }
    }

    using static_simplex_g3_type = Simplex<3, 3, 3>;
    auto const& staticSimplexTopology =
        geometry::GeometryPointTopology<
            static_simplex_g3_type>::get(
                RuntimeOrder{ 3 } );
    auto const& runtimeSimplexTopology =
        geometry::GeometryPointTopology<
            Simplex<3, Dynamic, 3>>::get(
                RuntimeOrder{ 3 } );
    for ( uint16_type face = 0;
          face < static_simplex_g3_type::numTopologicalFaces;
          ++face )
    {
        for ( uint16_type point = 0; point < 10; ++point )
        {
            BOOST_CHECK_EQUAL(
                staticSimplexTopology.faceToPoint(
                    face, point ),
                runtimeSimplexTopology.faceToPoint(
                    face, point ) );
            BOOST_CHECK_EQUAL(
                staticSimplexTopology.faceToPoint(
                    face, point ),
                static_simplex_g3_type::f2p(
                    face, point ) );
        }
    }

    using runtime_hypercube_mesh_type =
        Mesh<Hypercube<2, Dynamic>>;
    auto runtimeHypercubeMesh =
        runtime_hypercube_mesh_type::New( RuntimeOrder{ 4 } );
    auto dynamicHypercube =
        runtimeHypercubeMesh->makeElement();
    auto dynamicHypercubeFace =
        runtimeHypercubeMesh->makeFace();
    BOOST_CHECK_EQUAL( dynamicHypercube.order(), 4 );
    BOOST_CHECK_EQUAL( dynamicHypercube.nPoints(), 25 );
    BOOST_CHECK_EQUAL( dynamicHypercubeFace.order(), 4 );
    BOOST_CHECK_EQUAL( dynamicHypercubeFace.nPoints(), 5 );

    auto const& runtimeHexahedronTopology =
        geometry::GeometryPointTopology<
            Hypercube<3, Dynamic, 3>>::get(
                RuntimeOrder{ 4 } );
    BOOST_CHECK_EQUAL(
        runtimeHexahedronTopology.pointsOnFace(), 25 );
    BOOST_CHECK_EQUAL(
        runtimeHexahedronTopology.pointsOnEdge(), 5 );
    for ( uint16_type face = 0;
          face < Hypercube<3, Dynamic, 3>::numTopologicalFaces;
          ++face )
    {
        std::vector<bool> seen( 125, false );
        for ( uint16_type point = 0; point < 25; ++point )
        {
            auto const cellPoint =
                runtimeHexahedronTopology.faceToPoint(
                    face, point );
            BOOST_REQUIRE_LT( cellPoint, 125 );
            BOOST_CHECK( !seen[cellPoint] );
            seen[cellPoint] = true;
        }
    }

    GmshOrdering<Simplex<2, 2>> staticOrdering;
    GmshOrdering<Simplex<2, Dynamic>> runtimeOrdering{
        RuntimeOrder{ 2 } };
    BOOST_CHECK_EQUAL( runtimeOrdering.type(), staticOrdering.type() );
    BOOST_CHECK_EQUAL( runtimeOrdering.size(), staticOrdering.size() );
    for ( int point = 0; point < staticOrdering.size(); ++point )
    {
        BOOST_CHECK_EQUAL(
            runtimeOrdering.fromGmshId( point ),
            staticOrdering.fromGmshId( point ) );
        BOOST_CHECK_EQUAL(
            runtimeOrdering.toGmshId( point ),
            staticOrdering.toGmshId( point ) );
    }
    BOOST_CHECK_THROW(
        ( GmshOrdering<Simplex<2, Dynamic>>{
            RuntimeOrder{ 6 } } ),
        std::invalid_argument );

    GmshOrdering<Simplex<3, Dynamic>> tetraG5Ordering{
        RuntimeOrder{ 5 } };
    BOOST_CHECK_EQUAL( tetraG5Ordering.size(), 56 );
    BOOST_CHECK_EQUAL(
        tetraG5Ordering.fromGmshId( 4 ), 12 );
    BOOST_CHECK_EQUAL(
        tetraG5Ordering.fromGmshId( 28 ), 46 );
    BOOST_CHECK_EQUAL(
        tetraG5Ordering.fromGmshId( 51 ), 32 );

    GmshOrdering<Hypercube<2, Dynamic>>
        quadrangleG5Ordering{
            RuntimeOrder{ 5 } };
    BOOST_CHECK_EQUAL(
        quadrangleG5Ordering.type(),
        GMSH_QUADRANGLE_5 );
    BOOST_CHECK_EQUAL(
        quadrangleG5Ordering.size(), 36 );
    BOOST_CHECK_EQUAL(
        quadrangleG5Ordering.fromGmshId( 21 ),
        32 );
    BOOST_CHECK_EQUAL(
        quadrangleG5Ordering.fromGmshId( 35 ),
        26 );

    GmshOrdering<Hypercube<3, Dynamic>>
        hexahedronG3Ordering{
            RuntimeOrder{ 3 } };
    GmshOrdering<Hypercube<3, Dynamic>>
        hexahedronG5Ordering{
            RuntimeOrder{ 5 } };
    BOOST_CHECK_EQUAL(
        hexahedronG3Ordering.type(),
        GMSH_HEXAHEDRON_3 );
    BOOST_CHECK_EQUAL(
        hexahedronG3Ordering.size(), 64 );
    BOOST_CHECK_EQUAL(
        hexahedronG3Ordering.fromGmshId( 40 ),
        50 );
    BOOST_CHECK_EQUAL(
        hexahedronG5Ordering.type(),
        GMSH_HEXAHEDRON_5 );
    BOOST_CHECK_EQUAL(
        hexahedronG5Ordering.size(), 216 );
    BOOST_CHECK_EQUAL(
        hexahedronG5Ordering.fromGmshId( 153 ),
        200 );
    BOOST_CHECK_EQUAL(
        hexahedronG5Ordering.fromGmshId( 215 ),
        178 );

    std::stringstream archiveBuffer;
    {
        boost::archive::text_oarchive archive{
            archiveBuffer };
        archive << *runtimeMesh;
    }
    auto restoredRuntimeMesh =
        runtime_mesh_type::New( RuntimeOrder{ 1 } );
    {
        boost::archive::text_iarchive archive{
            archiveBuffer };
        archive >> *restoredRuntimeMesh;
    }
    BOOST_CHECK_EQUAL( restoredRuntimeMesh->order(), 6 );
    BOOST_CHECK(
        restoredRuntimeMesh->geometryRepresentation() ==
        geometry::GeometryRepresentationKind::BernsteinBezier );
    BOOST_REQUIRE( restoredRuntimeMesh->gm() );
    BOOST_CHECK_EQUAL(
        restoredRuntimeMesh->gm()->order(), 6 );
    BOOST_CHECK_EQUAL(
        restoredRuntimeMesh->geometryDescriptor()
            .coefficientCount(),
        28 );
}

/**
 * @test Compare static and runtime simplex mappings against exact G1--G3
 *       mappings in one, two, and three dimensions.
 */
BOOST_AUTO_TEST_CASE( analytic_simplex_static_runtime_g1_g2_g3 )
{
    checkAnalyticStaticRuntimeParity<1, 1, 1, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 1, 2, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 1, 3, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 2, 1, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 2, 2, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 2, 3, Simplex>();
    checkAnalyticStaticRuntimeParity<3, 3, 1, Simplex>();
    checkAnalyticStaticRuntimeParity<3, 3, 2, Simplex>();
    checkAnalyticStaticRuntimeParity<3, 3, 3, Simplex>();
}

/**
 * @test Compare static and runtime G1--G3 mappings for curves and surfaces
 *       embedded in a higher-dimensional physical space.
 */
BOOST_AUTO_TEST_CASE( analytic_embedded_static_runtime_g1_g2_g3 )
{
    checkAnalyticStaticRuntimeParity<1, 2, 1, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 2, 2, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 2, 3, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 3, 1, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 3, 2, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 3, 3, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 3, 1, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 3, 2, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 3, 3, Simplex>();
}

/**
 * @test Extend static/runtime analytic parity through the validated G4--G5
 *       envelope in 1D, 2D, 3D, embedded, and tensor-product settings.
 */
BOOST_AUTO_TEST_CASE( analytic_static_runtime_g4_g5 )
{
    checkAnalyticStaticRuntimeParity<1, 1, 4, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 1, 5, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 2, 4, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 2, 5, Simplex>();
    checkAnalyticStaticRuntimeParity<3, 3, 4, Simplex>();
    checkAnalyticStaticRuntimeParity<3, 3, 5, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 2, 4, Simplex>();
    checkAnalyticStaticRuntimeParity<1, 2, 5, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 3, 4, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 3, 5, Simplex>();
    checkAnalyticStaticRuntimeParity<2, 2, 4, Hypercube>();
    checkAnalyticStaticRuntimeParity<2, 2, 5, Hypercube>();
}

/**
 * @test Exercise true dynamic G6 maps above the configured static G1--G2
 *       mesh envelope in volume and embedded-surface settings.
 */
BOOST_AUTO_TEST_CASE( analytic_runtime_g6_above_static_envelope )
{
    checkAnalyticRuntimeOrder<2, 2, Simplex>( 6 );
    checkAnalyticRuntimeOrder<3, 3, Simplex>( 6 );
    checkAnalyticRuntimeOrder<2, 3, Simplex>( 6 );
    checkAnalyticRuntimeOrder<2, 2, Hypercube>( 6 );
}

/**
 * @test Verify curvature convergence on shrinking sphere, cylinder, and torus
 *       patches for G2 and G3, with static/runtime parity at every level.
 */
BOOST_AUTO_TEST_CASE( embedded_surface_curvature_refinement )
{
    auto checkOrder = []<int Order>()
    {
        std::array<double, 3> const sizes{ 0.4, 0.2, 0.1 };
        std::array<curvature_error_type, 3> sphereErrors;
        std::array<curvature_error_type, 3> cylinderErrors;
        std::array<curvature_error_type, 3> torusErrors;
        std::array<double, 3> sphereLaplaceErrors{};
        constexpr double majorRadius = 2.0;
        constexpr double minorRadius = 0.5;
        constexpr double barycentricCoordinate = -1.0/3.0;

        for ( std::size_t level = 0; level < sizes.size(); ++level )
        {
            double const h = sizes[level];
            sphereErrors[level] = surfaceCurvatureError<Order>(
                [h]( auto const& point )
                {
                    double const x = h*point( 0 );
                    double const y = h*point( 1 );
                    return Eigen::Vector3d(
                        x, y, std::sqrt( 1.0-x*x-y*y ) );
                },
                -1.0, 1.0 );

            withEmbeddedSurfaceContext<Order>(
                [h]( auto const& point )
                {
                    double const x = h*point( 0 );
                    double const y = h*point( 1 );
                    return Eigen::Vector3d(
                        x, y, std::sqrt( 1.0-x*x-y*y ) );
                },
                [&]( auto const& staticContext, auto const& runtimeContext )
                {
                    double const x = barycentricCoordinate;
                    double const y = barycentricCoordinate;
                    double const a = h*h;
                    double const z = std::sqrt( 1.0-a*( x*x+y*y ) );
                    Eigen::Vector2d gradient;
                    gradient << -a*x/z, -a*y/z;
                    Eigen::Matrix2d hessian;
                    hessian( 0, 0 ) = -a/z-a*a*x*x/( z*z*z );
                    hessian( 1, 1 ) = -a/z-a*a*y*y/( z*z*z );
                    hessian( 0, 1 ) = hessian( 1, 0 ) =
                        -a*a*x*y/( z*z*z );
                    double const staticLaplacian =
                        staticContext.laplaceBeltrami( 0, gradient, hessian );
                    double const runtimeLaplacian =
                        runtimeContext.laplaceBeltrami( 0, gradient, hessian );
                    BOOST_CHECK_SMALL( staticLaplacian-runtimeLaplacian,
                                       2e-10 );
                    sphereLaplaceErrors[level] =
                        std::abs( staticLaplacian+2.0*z );
                } );

            cylinderErrors[level] = surfaceCurvatureError<Order>(
                [h]( auto const& point )
                {
                    double const angle = h*point( 0 );
                    return Eigen::Vector3d(
                        std::cos( angle ), std::sin( angle ), h*point( 1 ) );
                },
                -0.5, 0.0 );

            double const minorAngle = h*barycentricCoordinate;
            double const firstPrincipal =
                -std::cos( minorAngle )/
                ( majorRadius+minorRadius*std::cos( minorAngle ) );
            double const secondPrincipal = -1.0/minorRadius;
            torusErrors[level] = surfaceCurvatureError<Order>(
                [h]( auto const& point )
                {
                    double const majorAngle = h*point( 0 );
                    double const minorAngleValue = h*point( 1 );
                    double const radial = majorRadius+
                                          minorRadius*std::cos( minorAngleValue );
                    return Eigen::Vector3d(
                        radial*std::cos( majorAngle ),
                        radial*std::sin( majorAngle ),
                        minorRadius*std::sin( minorAngleValue ) );
                },
                0.5*( firstPrincipal+secondPrincipal ),
                firstPrincipal*secondPrincipal );
        }

        auto checkReduction = []( auto const& errors )
        {
            for ( std::size_t level = 1; level < errors.size(); ++level )
            {
                BOOST_CHECK_LT( errors[level].mean,
                                0.8*errors[level-1].mean+1e-12 );
                BOOST_CHECK_LT( errors[level].gaussian,
                                0.8*errors[level-1].gaussian+1e-12 );
            }
        };
        checkReduction( sphereErrors );
        checkReduction( cylinderErrors );
        checkReduction( torusErrors );
        for ( std::size_t level = 1; level < sphereLaplaceErrors.size();
              ++level )
            BOOST_CHECK_LT( sphereLaplaceErrors[level],
                            0.8*sphereLaplaceErrors[level-1]+1e-12 );
    };

    checkOrder.template operator()<2>();
    checkOrder.template operator()<3>();
}

/**
 * @test Compare static and runtime G1--G3 mappings on a hypercube entity.
 */
BOOST_AUTO_TEST_CASE( analytic_hypercube_static_runtime_g1_g2_g3 )
{
    checkAnalyticStaticRuntimeParity<2, 2, 1, Hypercube>();
    checkAnalyticStaticRuntimeParity<2, 2, 2, Hypercube>();
    checkAnalyticStaticRuntimeParity<2, 2, 3, Hypercube>();
}

/**
 * @test Reject inverted and rank-deficient mappings and quantify a strongly
 *       anisotropic but nonsingular mapping.
 */
BOOST_AUTO_TEST_CASE( inverted_degenerate_and_nearly_singular_maps )
{
    withSimplex2DContext<1>(
        []( auto const& point )
        {
            return Eigen::Vector2d( -point( 0 ), point( 1 ) );
        },
        []( auto const& context )
        {
            BOOST_CHECK_LT( context.maximumSignedJacobian(), 0.0 );
            BOOST_CHECK( !context.isValid() );
            BOOST_CHECK( !context.mappingQuality( 0 ).orientation_preserving );
        } );

    withSimplex2DContext<1>(
        []( auto const& point )
        {
            return Eigen::Vector2d( point( 0 ), 0.0 );
        },
        []( auto const& context )
        {
            BOOST_CHECK_SMALL( context.minimumJacobian(), 1e-14 );
            BOOST_CHECK( !context.isValid() );
            BOOST_CHECK_SMALL( context.minimumSingularValue( 0 ), 1e-14 );
            BOOST_CHECK( std::isinf( context.conditionNumber( 0 ) ) );
        } );

    constexpr double epsilon = 1e-10;
    withSimplex2DContext<1>(
        []( auto const& point )
        {
            return Eigen::Vector2d( point( 0 ), epsilon * point( 1 ) );
        },
        []( auto const& context )
        {
            BOOST_CHECK( context.isValid() );
            BOOST_CHECK( !context.isValid( 1e-8 ) );
            BOOST_CHECK_GT( context.maximumConditionNumber(), 1e9 );
        } );
}

/**
 * @test Detect a curved mapping whose signed Jacobian changes sign and a
 *       mapping containing non-finite coordinates.
 */
BOOST_AUTO_TEST_CASE( curved_sign_change_and_nonfinite_detection )
{
    withSimplex2DContext<2>(
        []( auto const& point )
        {
            return Eigen::Vector2d( point( 0 ), point( 0 ) * point( 1 ) );
        },
        []( auto const& context )
        {
            BOOST_CHECK_LT( context.minimumSignedJacobian(), 0.0 );
            BOOST_CHECK_GT( context.maximumSignedJacobian(), 0.0 );
            BOOST_CHECK( !context.isValid() );
        } );

    withSimplex2DContext<1>(
        []( auto const& point )
        {
            return Eigen::Vector2d(
                point( 0 ) < 0 ? std::numeric_limits<double>::quiet_NaN() : point( 0 ),
                point( 1 ) );
        },
        []( auto const& context )
        {
            BOOST_CHECK( !context.mappingQuality( 0 ).finite );
            BOOST_CHECK( !context.isValid() );
        } );
}

BOOST_AUTO_TEST_SUITE_END()
