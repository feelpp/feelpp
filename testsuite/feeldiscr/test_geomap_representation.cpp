/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_geomap_representation.cpp
 * @brief Geometry representation, conversion, restriction, and subdivision tests.
 */

#define BOOST_TEST_MODULE geomap representation testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelpoly/concepts.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/geometryrepresentation.hpp>

#include <Eigen/Core>

#include <cmath>

using namespace Feel;

namespace
{
/** Type intentionally not satisfying the geometry representation contract. */
struct IncompleteRepresentation
{
};

/**
 * @brief Compare static/runtime nodal and Bezier representations on one entity.
 * @tparam Dim reference dimension
 * @tparam Order polynomial order
 * @tparam Entity reference-entity family
 */
template <int Dim, int Order, template <int, int, int> class Entity>
void checkRepresentationParity()
{
    using static_geomap_type = GeoMap<Dim, Order, Dim, double, Entity>;
    using runtime_geomap_type = GeoMap<Dim, Dynamic, Dim, double, Entity>;
    using static_nodal_type = geometry::NodalLagrange<Dim, Order, Entity, double>;
    using runtime_nodal_type = geometry::NodalLagrange<Dim, Dynamic, Entity, double>;
    using static_bezier_type = geometry::BernsteinBezier<Dim, Order, Entity, double>;
    using runtime_bezier_type = geometry::BernsteinBezier<Dim, Dynamic, Entity, double>;

    static_assert( GeometryRepresentation<static_nodal_type> );
    static_assert( GeometryRepresentation<runtime_nodal_type> );
    static_assert( GeometryRepresentation<static_bezier_type> );
    static_assert( GeometryRepresentation<runtime_bezier_type> );
    static_assert( !GeometryRepresentation<IncompleteRepresentation> );
    static_assert( NonLinearGeometricMappingConcept<
                   geometry::Mapping<Dim + 1, runtime_bezier_type>> );

    static_geomap_type staticGeomap;
    runtime_geomap_type runtimeGeomap{ RuntimeOrder{ Order } };
    auto staticNodal = staticGeomap.geometryRepresentation();
    auto runtimeNodal = runtimeGeomap.geometryRepresentation();
    static_bezier_type staticBezier;
    runtime_bezier_type runtimeBezier{ RuntimeOrder{ Order } };

    BOOST_REQUIRE_EQUAL( staticNodal.coefficientCount(), staticBezier.coefficientCount() );
    BOOST_REQUIRE_EQUAL( runtimeNodal.coefficientCount(), runtimeBezier.coefficientCount() );
    BOOST_CHECK( staticNodal.name() == "NodalLagrange" );
    BOOST_CHECK( staticBezier.name() == "BernsteinBezier" );
    BOOST_CHECK( static_geomap_type::coefficientSemantics() ==
                 geometry::CoefficientSemantics::InterpolationPoint );

    typename static_nodal_type::points_type evaluationPoints( Dim, 4 );
    for ( int d = 0; d < Dim; ++d )
        for ( Eigen::Index q = 0; q < evaluationPoints.cols(); ++q )
            evaluationPoints( d, q ) = -0.65 + 0.17 * d + 0.21 * q;
    if constexpr ( Feel::is_simplex_v<Entity<Dim, 1, Dim>> && Dim > 1 )
    {
        for ( Eigen::Index q = 0; q < evaluationPoints.cols(); ++q )
        {
            double sum = evaluationPoints.col( q ).sum();
            if ( sum > 2 - Dim )
                evaluationPoints.col( q ).array() -= ( sum - ( 2 - Dim ) ) / Dim;
        }
    }

    auto const staticNodalValues = staticNodal.basisValues( evaluationPoints );
    auto const runtimeNodalValues = runtimeNodal.basisValues( evaluationPoints );
    auto const staticBezierValues = staticBezier.basisValues( evaluationPoints );
    auto const runtimeBezierValues = runtimeBezier.basisValues( evaluationPoints );
    BOOST_CHECK_SMALL( ( staticNodalValues - runtimeNodalValues ).norm(), 2e-11 );
    BOOST_CHECK_SMALL( ( staticBezierValues - runtimeBezierValues ).norm(), 2e-13 );
    BOOST_CHECK_SMALL( ( staticNodal.basisValues( staticNodal.referenceNodes() ) -
                         Eigen::MatrixXd::Identity( staticNodal.coefficientCount(),
                                                    staticNodal.coefficientCount() ) )
                           .norm(),
                       2e-11 );
    for ( Eigen::Index q = 0; q < evaluationPoints.cols(); ++q )
    {
        BOOST_CHECK_SMALL( staticBezierValues.col( q ).sum() - 1.0, 2e-13 );
        typename static_geomap_type::node_t_type referencePoint( Dim );
        for ( int d = 0; d < Dim; ++d )
            referencePoint( d ) = evaluationPoints( d, q );
        for ( Eigen::Index basis = 0; basis < staticNodalValues.rows(); ++basis )
            BOOST_CHECK_SMALL( staticNodalValues( basis, q ) -
                                   staticGeomap.evaluate( basis, referencePoint )( 0 ),
                               2e-11 );
    }

    Eigen::MatrixXd nodalCoefficients( Dim + 1, staticNodal.coefficientCount() );
    for ( Eigen::Index component = 0; component < nodalCoefficients.rows(); ++component )
        for ( Eigen::Index coefficient = 0; coefficient < nodalCoefficients.cols(); ++coefficient )
            nodalCoefficients( component, coefficient ) =
                0.13 * ( component + 1 ) + 0.07 * coefficient + 0.005 * component * coefficient;
    auto const controlPoints = staticNodal.toBezier( nodalCoefficients );
    auto const recoveredNodal = staticNodal.fromBezier( controlPoints );
    BOOST_CHECK_SMALL( ( recoveredNodal - nodalCoefficients ).norm(), 2e-11 );
    BOOST_CHECK_SMALL( ( nodalCoefficients * staticNodalValues -
                         controlPoints * staticBezierValues )
                           .norm(),
                       2e-11 );

    auto const nodalGradients = staticNodal.basisGradients( evaluationPoints );
    auto const bezierGradients = staticBezier.basisGradients( evaluationPoints );
    for ( int d = 0; d < Dim; ++d )
    {
        BOOST_CHECK_SMALL( ( nodalCoefficients * nodalGradients[d] -
                             controlPoints * bezierGradients[d] )
                               .norm(),
                           5e-11 );
        for ( Eigen::Index q = 0; q < evaluationPoints.cols(); ++q )
            BOOST_CHECK_SMALL( bezierGradients[d].col( q ).sum(), 2e-12 );
    }

    auto const nodalHessians = staticNodal.basisHessians( evaluationPoints );
    auto const bezierHessians = staticBezier.basisHessians( evaluationPoints );
    for ( int derivative = 0; derivative < Dim * Dim; ++derivative )
    {
        BOOST_CHECK_SMALL( ( nodalCoefficients * nodalHessians[derivative] -
                             controlPoints * bezierHessians[derivative] )
                               .norm(),
                           2e-10 );
        for ( Eigen::Index q = 0; q < evaluationPoints.cols(); ++q )
            BOOST_CHECK_SMALL( bezierHessians[derivative].col( q ).sum(), 5e-12 );
    }

    geometry::CoefficientTable controlTable( staticBezier, controlPoints );
    static_assert( GeometryCoefficients<decltype( controlTable )> );
    BOOST_CHECK_EQUAL( controlTable.size(), staticBezier.coefficientCount() );
    geometry::CoefficientTable nodalTable( staticNodal, nodalCoefficients );
    auto convertedBezier = geometry::convert( nodalTable, staticBezier );
    auto convertedNodal = geometry::convert( convertedBezier, staticNodal );
    BOOST_CHECK_SMALL( ( convertedBezier.values() - controlPoints ).norm(), 2e-11 );
    BOOST_CHECK_SMALL( ( convertedNodal.values() - nodalCoefficients ).norm(), 2e-11 );

    geometry::Mapping<Dim + 1, static_nodal_type> nodalMapping( nodalTable );
    geometry::Mapping<Dim + 1, static_bezier_type> bezierMapping( convertedBezier );
    static_assert( GeometricMappingConcept<decltype( nodalMapping )> );
    static_assert( GeometricMappingConcept<decltype( bezierMapping )> );
    if constexpr ( Order == 1 )
    {
        static_assert( LinearGeometricMappingConcept<decltype( nodalMapping )> );
        static_assert( LinearGeometricMappingConcept<decltype( bezierMapping )> );
    }
    else
    {
        static_assert( NonLinearGeometricMappingConcept<decltype( nodalMapping )> );
        static_assert( NonLinearGeometricMappingConcept<decltype( bezierMapping )> );
    }
    BOOST_CHECK_SMALL( ( nodalMapping.transform( evaluationPoints ) -
                         bezierMapping.transform( evaluationPoints ) )
                           .norm(),
                       2e-11 );
    auto const nodalJacobians = nodalMapping.jacobians( evaluationPoints );
    auto const bezierJacobians = bezierMapping.jacobians( evaluationPoints );
    for ( size_type q = 0; q < nodalJacobians.size(); ++q )
        BOOST_CHECK_SMALL( ( nodalJacobians[q] - bezierJacobians[q] ).norm(), 5e-11 );
    auto const nodalMappingHessians = nodalMapping.hessians( evaluationPoints );
    auto const bezierMappingHessians = bezierMapping.hessians( evaluationPoints );
    for ( int derivative = 0; derivative < Dim * Dim; ++derivative )
        BOOST_CHECK_SMALL( ( nodalMappingHessians[derivative] -
                             bezierMappingHessians[derivative] )
                               .norm(),
                           2e-10 );
}
} // namespace

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( geomap_representation_suite )

/** @test Validate static/runtime representation and conversion parity. */
BOOST_AUTO_TEST_CASE( nodal_bezier_static_runtime_parity )
{
    checkRepresentationParity<1, 1, Simplex>();
    checkRepresentationParity<1, 3, Simplex>();
    checkRepresentationParity<2, 2, Simplex>();
    checkRepresentationParity<2, 3, Simplex>();
    checkRepresentationParity<2, 2, Hypercube>();
    checkRepresentationParity<2, 3, Hypercube>();
    checkRepresentationParity<3, 2, Simplex>();
    checkRepresentationParity<3, 3, Simplex>();
    checkRepresentationParity<3, 2, Hypercube>();
    checkRepresentationParity<3, 3, Hypercube>();
}

/** @test Check exact simplex and hypercube facet coefficient restrictions. */
BOOST_AUTO_TEST_CASE( exact_facet_restriction )
{
    constexpr int order = 3;
    using triangle_representation = geometry::BernsteinBezier<2, order, Simplex, double>;
    using segment_representation = geometry::BernsteinBezier<1, order, Simplex, double>;
    triangle_representation triangle;
    segment_representation segment;
    Eigen::MatrixXd triangleControls( 2, triangle.coefficientCount() );
    for ( Eigen::Index i = 0; i < triangleControls.cols(); ++i )
        triangleControls.col( i ) << 0.2 + 0.1 * i, -0.3 + 0.04 * i;
    geometry::CoefficientTable triangleTable( triangle, triangleControls );
    auto const facetControls = geometry::restrictToFacet(
        triangleTable, geometry::SimplexFacet{ 0 } );
    BOOST_REQUIRE_EQUAL( facetControls.cols(), order + 1 );

    typename triangle_representation::points_type trianglePoints( 2, 5 );
    typename segment_representation::points_type segmentPoints( 1, 5 );
    for ( Eigen::Index q = 0; q < 5; ++q )
    {
        double const parameter = -1.0 + 0.5 * q;
        trianglePoints.col( q ) << parameter, -parameter;
        segmentPoints( 0, q ) = -parameter;
    }
    BOOST_CHECK_SMALL( ( triangleControls * triangle.basisValues( trianglePoints ) -
                         facetControls * segment.basisValues( segmentPoints ) )
                           .norm(),
                       2e-12 );

    using square_representation = geometry::BernsteinBezier<2, order, Hypercube, double>;
    using edge_representation = geometry::BernsteinBezier<1, order, Hypercube, double>;
    square_representation square;
    edge_representation edge;
    Eigen::MatrixXd squareControls( 2, square.coefficientCount() );
    for ( Eigen::Index i = 0; i < squareControls.cols(); ++i )
        squareControls.col( i ) << -0.1 + 0.03 * i, 0.4 - 0.02 * i;
    geometry::CoefficientTable squareTable( square, squareControls );
    auto const edgeControls = geometry::restrictToFacet(
        squareTable, geometry::HypercubeFacet{ 0, true } );
    typename square_representation::points_type squarePoints( 2, 5 );
    typename edge_representation::points_type edgePoints( 1, 5 );
    for ( Eigen::Index q = 0; q < 5; ++q )
    {
        double const parameter = -1.0 + 0.5 * q;
        squarePoints.col( q ) << 1.0, parameter;
        edgePoints( 0, q ) = parameter;
    }
    BOOST_CHECK_SMALL( ( squareControls * square.basisValues( squarePoints ) -
                         edgeControls * edge.basisValues( edgePoints ) )
                           .norm(),
                       2e-12 );
}

/** @test Check de Casteljau subdivision against the original Bezier curve. */
BOOST_AUTO_TEST_CASE( bezier_curve_subdivision )
{
    using representation_type = geometry::BernsteinBezier<1, 3, Simplex, double>;
    representation_type representation;
    Eigen::MatrixXd controlPoints( 2, 4 );
    controlPoints << 0.0, 0.3, 0.8, 1.1,
        0.0, 0.7, -0.2, 0.4;
    geometry::CoefficientTable curve( representation, controlPoints );
    auto const [left, right] = geometry::subdivide( curve, 0.3 );

    representation_type::points_type childPoints( 1, 7 );
    representation_type::points_type leftParentPoints( 1, 7 );
    representation_type::points_type rightParentPoints( 1, 7 );
    for ( Eigen::Index q = 0; q < 7; ++q )
    {
        double const local = static_cast<double>( q ) / 6.0;
        childPoints( 0, q ) = 2.0 * local - 1.0;
        leftParentPoints( 0, q ) = 2.0 * ( 0.3 * local ) - 1.0;
        rightParentPoints( 0, q ) = 2.0 * ( 0.3 + 0.7 * local ) - 1.0;
    }
    BOOST_CHECK_SMALL( ( geometry::evaluate( left, childPoints ) -
                         geometry::evaluate( curve, leftParentPoints ) )
                           .norm(),
                       2e-12 );
    BOOST_CHECK_SMALL( ( geometry::evaluate( right, childPoints ) -
                         geometry::evaluate( curve, rightParentPoints ) )
                           .norm(),
                       2e-12 );
}

/** @test Reject invalid orders, coefficient cardinalities, and facet requests. */
BOOST_AUTO_TEST_CASE( representation_contract_diagnostics )
{
    using dynamic_bezier = geometry::BernsteinBezier<2, Dynamic, Simplex, double>;
    using static_bezier = geometry::BernsteinBezier<2, 2, Simplex, double>;
    BOOST_CHECK_THROW( dynamic_bezier( RuntimeOrder{ 0 } ), std::invalid_argument );
    BOOST_CHECK_THROW( static_bezier( RuntimeOrder{ 3 } ), std::invalid_argument );

    static_bezier bezier;
    Eigen::MatrixXd wrongCardinality( 2, bezier.coefficientCount() - 1 );
    BOOST_CHECK_THROW( geometry::CoefficientTable( bezier, wrongCardinality ),
                       std::invalid_argument );
    BOOST_CHECK_THROW( static_cast<void>( bezier.facetCoefficientIndices(
                           geometry::SimplexFacet{ 3 } ) ),
                       std::out_of_range );

    Eigen::MatrixXd controls( 2, bezier.coefficientCount() );
    controls.setZero();
    geometry::CoefficientTable coefficientTable( bezier, controls );
    BOOST_CHECK_THROW( ( geometry::Mapping<3, static_bezier>( coefficientTable ) ),
                       std::invalid_argument );

    using curve_representation = geometry::BernsteinBezier<1, 2, Simplex, double>;
    curve_representation curveRepresentation;
    Eigen::MatrixXd curveControls( 2, curveRepresentation.coefficientCount() );
    curveControls.setZero();
    geometry::CoefficientTable curve( curveRepresentation, curveControls );
    BOOST_CHECK_THROW( static_cast<void>( geometry::subdivide( curve, -0.1 ) ),
                       std::out_of_range );
    BOOST_CHECK_THROW( static_cast<void>( geometry::subdivide( curve, 1.1 ) ),
                       std::out_of_range );
}

BOOST_AUTO_TEST_SUITE_END()
