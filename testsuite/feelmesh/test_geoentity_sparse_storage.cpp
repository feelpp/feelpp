/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_geoentity_sparse_storage.cpp
 * @brief Validate compact geometric-entity cores and mesh-owned sidecars.
 */

#define BOOST_TEST_MAIN
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/points.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <array>
#include <limits>
#include <sstream>
#include <type_traits>

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
using entity_type = Feel::GeoEntity<Feel::Simplex<2, 1>>;
using core_type = typename entity_type::core_type;
using storage_type = typename entity_type::sparse_metadata_storage_type;

static_assert( std::is_standard_layout_v<core_type> );
static_assert( std::is_trivially_copyable_v<core_type> );
static_assert( sizeof( core_type ) <= 32 );
} // namespace

/**
 * @brief Check independent sparse payload allocation and entity rekeying.
 */
BOOST_AUTO_TEST_CASE( compact_core_and_independent_sidecars )
{
    entity_type entity( 42 );
    storage_type storage;

    entity.attachSparseMetadataStorage( storage );
    BOOST_CHECK( entity.hasSparseMetadataStorage() );
    BOOST_CHECK_EQUAL( storage.statistics().payloadCount(), 0 );

    entity.setMarker( 1, 17 );
    auto statistics = storage.statistics();
    BOOST_CHECK_EQUAL( statistics.markerEntries, 1 );
    BOOST_CHECK_EQUAL( statistics.remoteIdEntries, 0 );
    BOOST_CHECK_EQUAL( statistics.incidenceEntries, 0 );
    BOOST_CHECK_EQUAL( statistics.measureEntries, 0 );

    entity.setIdInOtherPartitions( 3, 91 );
    entity.addElement( 12, 2 );
    entity.mutableMeasureMetadata().measure = 2.5;
    statistics = storage.statistics();
    BOOST_CHECK_EQUAL( statistics.payloadCount(), 4 );
    BOOST_CHECK_EQUAL( entity.idInOthersPartitions( 3 ), 91 );
    BOOST_CHECK_EQUAL( entity.numberOfElements(), 1 );
    BOOST_CHECK_CLOSE( *entity.measureMetadata().measure, 2.5, 1e-12 );

    std::array<entity_type::index_type, 3> connectivity{
        4, 8, 15
    };
    entity.setConnectivity( connectivity );
    BOOST_CHECK_NE(
        entity.connectivityOffset(),
        std::numeric_limits<Feel::uint32_type>::max() );
    BOOST_CHECK_EQUAL( entity.connectivityCount(), 3 );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        entity.connectivityIds().begin(),
        entity.connectivityIds().end(),
        connectivity.begin(), connectivity.end() );

    entity.setId( 43 );
    BOOST_CHECK_EQUAL( entity.id(), 43 );
    BOOST_CHECK_EQUAL( entity.idInOthersPartitions( 3 ), 91 );
    BOOST_CHECK_EQUAL( storage.statistics().payloadCount(), 4 );
}

/**
 * @brief Reject invalid and repeated compact-connectivity releases.
 */
BOOST_AUTO_TEST_CASE( connectivity_storage_checks_releases )
{
    storage_type storage;
    std::array<entity_type::index_type, 3> connectivity{
        4, 8, 15
    };
    auto const offset =
        storage.allocateConnectivity( connectivity );

    BOOST_CHECK_THROW(
        storage.releaseConnectivity( offset + 2, 3 ),
        std::out_of_range );
    storage.releaseConnectivity( offset, 3 );
    BOOST_CHECK_THROW(
        storage.releaseConnectivity( offset, 3 ),
        std::logic_error );
    auto const statistics =
        storage.connectivityStatistics();
    BOOST_CHECK_EQUAL( statistics.activeEntries, 0 );
    BOOST_CHECK_EQUAL( statistics.reusableEntries, 3 );
}

/**
 * @brief Check that copies do not alias another mesh's sparse metadata.
 */
BOOST_AUTO_TEST_CASE( copy_and_move_preserve_metadata_ownership )
{
    entity_type entity( 7 );
    storage_type storage;
    entity.attachSparseMetadataStorage( storage );
    entity.setMarker( 1, 5 );

    entity_type copy( entity );
    BOOST_CHECK( !copy.hasSparseMetadataStorage() );
    BOOST_CHECK_EQUAL( copy.marker( 1 ).value(), 5 );
    copy.setMarker( 1, 9 );
    BOOST_CHECK_EQUAL( entity.marker( 1 ).value(), 5 );
    BOOST_CHECK_EQUAL( copy.marker( 1 ).value(), 9 );

    entity_type moved( std::move( entity ) );
    BOOST_CHECK( moved.hasSparseMetadataStorage() );
    BOOST_CHECK( !entity.hasSparseMetadataStorage() );
    BOOST_CHECK_EQUAL( moved.marker( 1 ).value(), 5 );

    storage_type destinationStorage;
    moved.attachSparseMetadataStorage( destinationStorage );
    BOOST_CHECK_EQUAL( storage.statistics().payloadCount(), 0 );
    BOOST_CHECK_EQUAL(
        destinationStorage.statistics().markerEntries, 1 );
    BOOST_CHECK_EQUAL( moved.marker( 1 ).value(), 5 );
}

/**
 * @brief Check point-container transfer, copy isolation, erase, and restart.
 */
BOOST_AUTO_TEST_CASE( point_container_owns_sparse_metadata )
{
    using points_type = Feel::Points<2>;
    using point_type = typename points_type::point_type;

    points_type points;
    point_type point( 10 );
    point( 0 ) = 1.25;
    point( 1 ) = -2.5;
    point.setMarker( 1, 23 );
    point.setIdInOtherPartitions( 2, 101 );
    point.addElement( 4, 0 );

    auto [iterator, inserted] = points.addPoint( std::move( point ) );
    BOOST_REQUIRE( inserted );
    BOOST_CHECK( iterator->second.hasSparseMetadataStorage() );
    auto statistics = points.pointSparseMetadataStatistics();
    BOOST_CHECK_EQUAL( statistics.markerEntries, 1 );
    BOOST_CHECK_EQUAL( statistics.remoteIdEntries, 1 );
    BOOST_CHECK_EQUAL( statistics.incidenceEntries, 1 );
    BOOST_CHECK_EQUAL( statistics.measureEntries, 0 );
    auto coordinates = points.pointCoordinates();
    BOOST_CHECK_EQUAL( coordinates.rows(), 2 );
    BOOST_CHECK_EQUAL( coordinates.cols(), 1 );
    BOOST_CHECK_CLOSE( coordinates( 0, 0 ), 1.25, 1e-12 );
    BOOST_CHECK_CLOSE( coordinates( 1, 0 ), -2.5, 1e-12 );
    BOOST_CHECK_EQUAL( points.pointLocalIndex( 10 ), 0 );
    BOOST_CHECK_EQUAL(
        points.pointCoordinateStorageStatistics()
            .externalLookupEntries,
        0 );

    points.editPointCoordinates(
        []( auto& coordinateMap )
        {
            coordinateMap.col( 0 ) *= 2.;
        } );
    BOOST_CHECK_CLOSE( points.point( 10 )( 0 ), 2.5, 1e-12 );
    BOOST_CHECK_CLOSE( points.point( 10 )( 1 ), -5., 1e-12 );

    points_type copied( points );
    copied.point( 10 ).setMarker( 1, 29 );
    BOOST_CHECK_EQUAL( points.point( 10 ).marker( 1 ).value(), 23 );
    BOOST_CHECK_EQUAL( copied.point( 10 ).marker( 1 ).value(), 29 );

    std::stringstream stream(
        std::ios::in | std::ios::out | std::ios::binary );
    {
        boost::archive::binary_oarchive archive( stream );
        archive << points;
    }

    points_type restored;
    {
        boost::archive::binary_iarchive archive( stream );
        archive >> restored;
    }
    BOOST_REQUIRE( restored.hasPoint( 10 ) );
    BOOST_CHECK( restored.point( 10 ).hasSparseMetadataStorage() );
    BOOST_CHECK_EQUAL( restored.point( 10 ).marker( 1 ).value(), 23 );
    BOOST_CHECK_EQUAL(
        restored.point( 10 ).idInOthersPartitions( 2 ), 101 );

    points.erasePoint( points.pointIterator( 10 ) );
    BOOST_CHECK_EQUAL(
        points.pointSparseMetadataStatistics().payloadCount(), 0 );
    BOOST_CHECK_EQUAL(
        points.pointCoordinateStorageStatistics().points, 0 );
}

/**
 * @brief Check that dense point IDs require no external lookup table.
 */
BOOST_AUTO_TEST_CASE( dense_point_coordinates_are_eigen_mapped )
{
    using points_type = Feel::Points<3>;
    using point_type = typename points_type::point_type;

    points_type points;
    point_type point0( 0 );
    point_type point1( 1 );
    point0( 0 ) = 1.;
    point0( 1 ) = 2.;
    point0( 2 ) = 3.;
    point1( 0 ) = 4.;
    point1( 1 ) = 5.;
    point1( 2 ) = 6.;
    BOOST_REQUIRE( points.addPoint( std::move( point0 ) ).second );
    BOOST_REQUIRE( points.addPoint( std::move( point1 ) ).second );

    auto const statistics =
        points.pointCoordinateStorageStatistics();
    BOOST_CHECK_EQUAL( statistics.points, 2 );
    BOOST_CHECK_EQUAL( statistics.scalarEntries, 6 );
    BOOST_CHECK_EQUAL( statistics.externalLookupEntries, 0 );

    auto coordinates =
        static_cast<points_type const&>( points ).pointCoordinates();
    BOOST_CHECK_EQUAL( coordinates.rows(), 3 );
    BOOST_CHECK_EQUAL( coordinates.cols(), 2 );
    BOOST_CHECK_CLOSE( coordinates( 2, 1 ), 6., 1e-12 );

    point_type point3( 3 );
    point3( 0 ) = 7.;
    point3( 1 ) = 8.;
    point3( 2 ) = 9.;
    BOOST_REQUIRE( points.addPoint( std::move( point3 ) ).second );
    BOOST_CHECK_EQUAL( points.pointLocalIndex( 3 ), 2 );
    BOOST_CHECK_EQUAL(
        points.pointCoordinateStorageStatistics()
            .externalLookupEntries,
        3 );
}

/**
 * @brief Check that cached cell measures occupy only the measure sidecar.
 */
BOOST_AUTO_TEST_CASE( element_measure_cache_is_sparse )
{
    using elements_type = Feel::Elements<Feel::Simplex<2, 1>>;
    using element_type = typename elements_type::element_type;
    using point_type = typename element_type::point_type;

    elements_type elements;
    element_type element( 31 );
    point_type point0( 10 );
    point_type point1( 11 );
    point_type point2( 12 );
    element.setPoint( 0, point0 );
    element.setPoint( 1, point1 );
    element.setPoint( 2, point2 );
    element.setMeasure( 3.25 );
    element.setFaceMeasure( 0, 1.5 );

    auto [iterator, inserted] =
        elements.addElement( std::move( element ) );
    BOOST_REQUIRE( inserted );
    BOOST_CHECK( iterator->second.hasSparseMetadataStorage() );
    auto const statistics =
        elements.elementSparseMetadataStatistics();
    BOOST_CHECK_EQUAL( statistics.measureEntries, 1 );
    BOOST_CHECK_EQUAL( statistics.markerEntries, 0 );
    BOOST_CHECK_EQUAL( statistics.incidenceEntries, 0 );
    BOOST_CHECK_EQUAL( statistics.remoteIdEntries, 0 );
    BOOST_CHECK_CLOSE( iterator->second.measure(), 3.25, 1e-12 );
    BOOST_CHECK_CLOSE( iterator->second.faceMeasure( 0 ), 1.5, 1e-12 );
    std::array<element_type::index_type, 3> expectedIds{
        10, 11, 12
    };
    BOOST_CHECK_EQUAL_COLLECTIONS(
        iterator->second.pointIds().begin(),
        iterator->second.pointIds().end(),
        expectedIds.begin(), expectedIds.end() );
    auto const connectivityStatistics =
        elements.elementConnectivityStatistics();
    BOOST_CHECK_EQUAL( connectivityStatistics.allocatedEntries, 3 );
    BOOST_CHECK_EQUAL( connectivityStatistics.activeEntries, 3 );
    BOOST_CHECK_EQUAL( connectivityStatistics.reusableEntries, 0 );

    auto const firstOffset = iterator->second.connectivityOffset();
    elements.eraseElementOnly( iterator );
    BOOST_CHECK_EQUAL(
        elements.elementConnectivityStatistics().activeEntries, 0 );
    BOOST_CHECK_EQUAL(
        elements.elementConnectivityStatistics().reusableEntries, 3 );

    element_type replacement( 32 );
    replacement.setPoint( 0, point2 );
    replacement.setPoint( 1, point1 );
    replacement.setPoint( 2, point0 );
    auto [replacementIterator, replacementInserted] =
        elements.addElement( std::move( replacement ) );
    BOOST_REQUIRE( replacementInserted );
    BOOST_CHECK_EQUAL(
        replacementIterator->second.connectivityOffset(),
        firstOffset );
}

/**
 * @brief Check compact connectivity and coordinates on a production mesh.
 */
BOOST_AUTO_TEST_CASE( compact_mesh_geometry_matches_legacy_facade )
{
    using mesh_type = Feel::Mesh<Feel::Simplex<2, 1>>;

    auto mesh = Feel::createGMSHMesh(
        Feel::_mesh = new mesh_type,
        Feel::_desc = Feel::domain(
            Feel::_name = "compact-mesh-geometry",
            Feel::_shape = "hypercube", Feel::_dim = 2,
            Feel::_h = 0.4 ),
        Feel::_verbose = 0 );

    Feel::size_type expectedElementEntries = 0;
    for ( auto const& [id, element] : mesh->elements() )
    {
        (void)id;
        auto const pointIds = element.pointIds();
        BOOST_REQUIRE_EQUAL(
            pointIds.size(), element.nPoints() );
        expectedElementEntries += pointIds.size();
        for ( Feel::uint16_type localPoint = 0;
              localPoint < pointIds.size(); ++localPoint )
            BOOST_CHECK_EQUAL(
                pointIds[localPoint],
                element.point( localPoint ).id() );
    }
    BOOST_CHECK_EQUAL(
        mesh->elementConnectivityStatistics().activeEntries,
        expectedElementEntries );

    Feel::size_type expectedFaceEntries = 0;
    for ( auto const& [id, face] : mesh->faces() )
    {
        (void)id;
        auto const pointIds = face.pointIds();
        expectedFaceEntries += pointIds.size();
        for ( Feel::uint16_type localPoint = 0;
              localPoint < pointIds.size(); ++localPoint )
            BOOST_CHECK_EQUAL(
                pointIds[localPoint],
                face.point( localPoint ).id() );
    }
    BOOST_CHECK_EQUAL(
        mesh->faceConnectivityStatistics().activeEntries,
        expectedFaceEntries );

    auto coordinates =
        static_cast<mesh_type const&>( *mesh ).pointCoordinates();
    BOOST_REQUIRE_EQUAL( coordinates.cols(), mesh->numPoints() );
    for ( auto const& [id, point] : mesh->points() )
    {
        auto const localIndex = mesh->pointLocalIndex( id );
        for ( Feel::uint16_type component = 0;
              component < mesh_type::nRealDim; ++component )
            BOOST_CHECK_CLOSE(
                coordinates( component, localIndex ),
                point( component ), 1e-12 );
    }
    BOOST_CHECK_EQUAL(
        mesh->pointCoordinateStorageStatistics()
            .externalLookupEntries,
        0 );

    auto pointIterator = mesh->beginPoint();
    BOOST_REQUIRE( pointIterator != mesh->endPoint() );
    auto const movedPointId = pointIterator->first;
    pointIterator->second( 0 ) += 0.125;
    auto const movedCoordinate = pointIterator->second( 0 );
    mesh->updateForUseAfterMovingNodes( false );
    auto updatedCoordinates =
        static_cast<mesh_type const&>( *mesh ).pointCoordinates();
    BOOST_CHECK_CLOSE(
        updatedCoordinates(
            0, mesh->pointLocalIndex( movedPointId ) ),
        movedCoordinate, 1e-12 );
}
