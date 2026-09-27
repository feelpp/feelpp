/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2012-2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE mesh_face_cleanup
#include <feel/feelcore/testsuite.hpp>
#include <feel/feeldiscr/mesh.hpp>

#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

using namespace Feel;
using SurfaceMesh = Mesh<Simplex<2, 1, 3>>;
using VolumeMesh = Mesh<Simplex<3>>;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** @brief Make a sequential surface mesh for face cleanup checks. */
std::shared_ptr<SurfaceMesh> makeMesh()
{
    return std::make_shared<SurfaceMesh>( Environment::worldCommSeqPtr() );
}

/** @brief Add a face with a distinct marker to verify retained identity. */
void addFace( SurfaceMesh& mesh, SurfaceMesh::size_type id )
{
    SurfaceMesh::face_type face;
    face.setId( id );
    face.setMarker( id + 100 );
    mesh.addFace( std::move( face ) );
}

/** @brief Read face IDs in their stored ordered sequence. */
std::vector<SurfaceMesh::size_type> orderedIds( SurfaceMesh const& mesh )
{
    std::vector<SurfaceMesh::size_type> ids;
    for ( auto const& face : mesh.orderedFaces() )
        ids.push_back( unwrap_ref( face ).id() );
    return ids;
}

/** @brief Add a two-dimensional point to the surface test mesh. */
void addPoint( SurfaceMesh& mesh, SurfaceMesh::size_type id, double x, double y )
{
    SurfaceMesh::point_type point( id );
    point( 0 ) = x;
    point( 1 ) = y;
    point( 2 ) = 0;
    point.setProcessId( 0 );
    point.setProcessIdInPartition( 0 );
    mesh.addPoint( std::move( point ) );
}

/** @brief Check that face cleanup preserves triangle connectivity and markers. */
void checkConnectivityCleanup( size_type updateFlags )
{
    auto mesh = makeMesh();
    addPoint( *mesh, 0, 0, 0 );
    addPoint( *mesh, 1, 1, 0 );
    addPoint( *mesh, 2, 0, 1 );
    addPoint( *mesh, 3, 1, 1 );
    std::array<std::array<SurfaceMesh::size_type, 3>, 2> triangles = {{{0, 1, 2}, {1, 3, 2}}};
    for ( SurfaceMesh::size_type id = 0; id < triangles.size(); ++id )
    {
        SurfaceMesh::element_type element;
        element.setId( id );
        element.setMarker( 50 );
        element.setProcessId( 0 );
        element.setProcessIdInPartition( 0 );
        for ( uint16_type vertex = 0; vertex < 3; ++vertex )
            element.setPoint( vertex, mesh->point( triangles[id][vertex] ) );
        mesh->addElement( std::move( element ) );
    }
    std::array<std::array<SurfaceMesh::size_type, 2>, 5> edges = {{{0, 1}, {1, 2}, {2, 0}, {1, 3}, {3, 2}}};
    std::map<SurfaceMesh::size_type, SurfaceMesh::face_type const*> addresses;
    for ( SurfaceMesh::size_type id = 0; id < edges.size(); ++id )
    {
        SurfaceMesh::face_type face;
        face.setId( 100 + id );
        face.setMarker( 11 + id );
        face.setProcessId( 0 );
        face.setProcessIdInPartition( 0 );
        for ( uint16_type vertex = 0; vertex < 2; ++vertex )
            face.setPoint( vertex, mesh->point( edges[id][vertex] ) );
        auto inserted = mesh->addFace( std::move( face ) );
        addresses.emplace( 100 + id, &inserted.first->second );
    }
    // Distinct disconnected node pairs reproduce the Chicago cleanup, without
    // introducing duplicate edges or removing any surface triangles.
    for ( SurfaceMesh::size_type id = 0; id < 4096; ++id )
    {
        auto firstPoint = 10 + 2 * id;
        addPoint( *mesh, firstPoint, 10 + id, 0 );
        addPoint( *mesh, firstPoint + 1, 10 + id, 1 );
        SurfaceMesh::face_type face;
        face.setId( 1000 + id );
        face.setMarker( 99 );
        face.setProcessId( 0 );
        face.setProcessIdInPartition( 0 );
        face.setPoint( 0, mesh->point( firstPoint ) );
        face.setPoint( 1, mesh->point( firstPoint + 1 ) );
        mesh->addFace( std::move( face ) );
        if ( updateFlags & MESH_UPDATE_FACES_MINIMAL )
        {
            SurfaceMesh::face_type duplicate;
            duplicate.setId( 10000 + id );
            duplicate.setMarker( 98 );
            duplicate.setProcessId( 0 );
            duplicate.setProcessIdInPartition( 0 );
            duplicate.setPoint( 0, mesh->point( firstPoint ) );
            duplicate.setPoint( 1, mesh->point( firstPoint + 1 ) );
            mesh->addFace( std::move( duplicate ) );
        }
    }
    mesh->components().reset();
    mesh->components().set( updateFlags | MESH_NO_UPDATE_MEASURES );
    mesh->updateForUse();
    BOOST_CHECK_EQUAL( mesh->elements().size(), 2 );
    BOOST_REQUIRE_EQUAL( mesh->faces().size(), 5 );
    BOOST_REQUIRE_EQUAL( mesh->orderedFaces().size(), 5 );
    for ( auto const& [id, address] : addresses )
    {
        auto const& face = mesh->faces().at( id );
        BOOST_CHECK( &face == address );
        BOOST_CHECK_EQUAL( face.marker().value(), 11 + id - 100 );
        BOOST_CHECK( face.isConnectedTo0() );
        BOOST_CHECK_EQUAL( face.isConnectedTo1(), id == 101 );
        BOOST_CHECK_EQUAL( face.isOnBoundary(), id != 101 );
    }
    for ( auto const& [id, element] : mesh->elements() )
    {
        BOOST_CHECK_EQUAL( element.marker().value(), 50 );
        for ( uint16_type vertex = 0; vertex < 3; ++vertex )
            BOOST_CHECK_EQUAL( element.point( vertex ).id(), triangles[id][vertex] );
    }
}

/** @brief Check that a registered duplicate face retains the first face. */
void checkConnectedFaceDuplicate( size_type updateFlags )
{
    auto mesh = makeMesh();
    addPoint( *mesh, 0, 0, 0 );
    addPoint( *mesh, 1, 1, 0 );
    addPoint( *mesh, 2, 0, 1 );
    SurfaceMesh::element_type element;
    element.setId( 0 );
    element.setProcessId( 0 );
    element.setProcessIdInPartition( 0 );
    for ( uint16_type vertex = 0; vertex < 3; ++vertex )
        element.setPoint( vertex, mesh->point( vertex ) );
    mesh->addElement( std::move( element ) );

    for ( SurfaceMesh::size_type id : {100, 101} )
    {
        SurfaceMesh::face_type face;
        face.setId( id );
        face.setMarker( id + 10 );
        face.setProcessId( 0 );
        face.setProcessIdInPartition( 0 );
        face.setPoint( 0, mesh->point( 0 ) );
        face.setPoint( 1, mesh->point( 1 ) );
        mesh->addFace( std::move( face ) );
    }
    auto const* retained = &mesh->beginFace()->second;
    auto const retainedId = retained->id();
    auto const retainedMarker = retained->marker().value();

    mesh->components().reset();
    mesh->components().set( updateFlags | MESH_NO_UPDATE_MEASURES );
    mesh->updateForUse();

    BOOST_REQUIRE_EQUAL( mesh->faces().size(), 3 );
    BOOST_CHECK( &mesh->faces().at( retainedId ) == retained );
    BOOST_CHECK_EQUAL( retained->marker().value(), retainedMarker );
    BOOST_CHECK( retained->isConnectedTo0() );
    bool elementUsesRetainedFace = false;
    auto const& updatedElement = mesh->elements().at( 0 );
    for ( uint16_type face = 0; face < updatedElement.numLocalEdges; ++face )
        elementUsesRetainedFace |= updatedElement.facePtr( face ) == retained;
    BOOST_CHECK( elementUsesRetainedFace );
}

/** @brief Add a three-dimensional point to the volume test mesh. */
void addVolumePoint( VolumeMesh& mesh, VolumeMesh::size_type id, double x, double y, double z )
{
    VolumeMesh::point_type point( id );
    point( 0 ) = x;
    point( 1 ) = y;
    point( 2 ) = z;
    point.setProcessId( 0 );
    point.setProcessIdInPartition( 0 );
    mesh.addPoint( std::move( point ) );
}
}

/** @brief Bulk face erasure retains surviving addresses and ordered IDs. */
BOOST_AUTO_TEST_CASE( bulkErasePreservesSurvivors )
{
    auto mesh = makeMesh();
    for ( auto id : {40u, 3u, 20u, 1u, 9u, 60u} )
        addFace( *mesh, id );
    auto const* retained = &mesh->faces().at( 20 );
    std::size_t calls = 0;
    auto const removed = mesh->eraseFacesIf( [&calls]( auto const& face ) { ++calls; return face.id() % 3 == 0; } );
    BOOST_CHECK_EQUAL( removed, 3 );
    BOOST_CHECK_EQUAL( calls, 6 );
    BOOST_CHECK( &mesh->faces().at( 20 ) == retained );
    BOOST_CHECK_EQUAL( retained->marker().value(), 120 );
    BOOST_CHECK( orderedIds( *mesh ) == ( std::vector<SurfaceMesh::size_type>{40, 20, 1} ) );
    mesh->updateOrderedFaces();
    BOOST_CHECK( orderedIds( *mesh ) == ( std::vector<SurfaceMesh::size_type>{1, 20, 40} ) );
    BOOST_CHECK_EQUAL( mesh->eraseFacesIf( []( auto const& ) { return false; } ), 0 );
    BOOST_CHECK_EQUAL( mesh->eraseFacesIf( []( auto const& ) { return true; } ), 3 );
    BOOST_CHECK( mesh->faces().empty() );
    BOOST_CHECK( mesh->orderedFaces().empty() );
    BOOST_CHECK_EQUAL( mesh->eraseFacesIf( []( auto const& ) { return true; } ), 0 );
}

/** @brief A throwing erase predicate leaves the mesh unchanged. */
BOOST_AUTO_TEST_CASE( throwingPredicateDoesNotMutateMesh )
{
    auto mesh = makeMesh();
    for ( auto id : {7u, 3u, 9u} )
        addFace( *mesh, id );
    auto const before = orderedIds( *mesh );
    auto const* retained = &mesh->faces().at( 7 );
    std::size_t calls = 0;
    auto predicate = [&calls]( auto const& )
    {
        if ( ++calls == 2 )
            throw std::runtime_error( "predicate failure" );
        return true;
    };
    BOOST_CHECK_THROW( mesh->eraseFacesIf( predicate ), std::runtime_error );
    BOOST_CHECK_EQUAL( mesh->faces().size(), 3 );
    BOOST_CHECK( orderedIds( *mesh ) == before );
    BOOST_CHECK( &mesh->faces().at( 7 ) == retained );
}

/** @brief Minimal face updates keep element topology and registered markers. */
BOOST_AUTO_TEST_CASE( minimalConnectivityRetainsMarkersAndTriangles )
{
    checkConnectivityCleanup( MESH_UPDATE_FACES_MINIMAL );
}

/** @brief Full face updates keep element topology and registered markers. */
BOOST_AUTO_TEST_CASE( fullConnectivityRetainsMarkersAndTriangles )
{
    checkConnectivityCleanup( MESH_UPDATE_FACES );
}

/** @brief Minimal updates preserve the first registered duplicate face. */
BOOST_AUTO_TEST_CASE( connectedDuplicateFaceMinimal )
{
    checkConnectedFaceDuplicate( MESH_UPDATE_FACES_MINIMAL );
}

/** @brief Full updates preserve the first registered duplicate face. */
BOOST_AUTO_TEST_CASE( connectedDuplicateFaceFull )
{
    checkConnectedFaceDuplicate( MESH_UPDATE_FACES );
}

/** @brief Volume edge cleanup preserves the first edge and element links. */
BOOST_AUTO_TEST_CASE( connectedDuplicateEdgesPreserveFirstEdge )
{
    auto mesh = std::make_shared<VolumeMesh>( Environment::worldCommSeqPtr() );
    addVolumePoint( *mesh, 0, 0, 0, 0 );
    addVolumePoint( *mesh, 1, 1, 0, 0 );
    addVolumePoint( *mesh, 2, 0, 1, 0 );
    addVolumePoint( *mesh, 3, 0, 0, 1 );
    VolumeMesh::element_type element;
    element.setId( 0 );
    element.setProcessId( 0 );
    element.setProcessIdInPartition( 0 );
    for ( uint16_type vertex = 0; vertex < 4; ++vertex )
        element.setPoint( vertex, mesh->point( vertex ) );
    mesh->addElement( std::move( element ) );

    for ( VolumeMesh::size_type index = 0; index < 4096; ++index )
    {
        VolumeMesh::edge_type edge;
        edge.setId( 100 + index );
        edge.setMarker( index + 10 );
        edge.setProcessId( 0 );
        edge.setProcessIdInPartition( 0 );
        edge.setPoint( 0, mesh->point( 0 ) );
        edge.setPoint( 1, mesh->point( 1 ) );
        mesh->addEdge( std::move( edge ) );
    }
    auto const* retained = &mesh->beginEdge()->second;
    auto const retainedId = retained->id();
    auto const retainedMarker = retained->marker().value();

    mesh->components().reset();
    mesh->components().set( MESH_UPDATE_FACES_MINIMAL | MESH_UPDATE_EDGES | MESH_NO_UPDATE_MEASURES );
    mesh->updateForUse();

    BOOST_REQUIRE_EQUAL( mesh->edges().size(), 6 );
    BOOST_REQUIRE_EQUAL( std::distance( mesh->beginOrderedEdge(), mesh->endOrderedEdge() ), 6 );
    BOOST_CHECK( &mesh->edges().at( retainedId ) == retained );
    BOOST_CHECK_EQUAL( retained->marker().value(), retainedMarker );
    auto const& updatedElement = mesh->elements().at( 0 );
    bool elementUsesRetainedEdge = false;
    for ( uint16_type edge = 0; edge < updatedElement.numEdges; ++edge )
        elementUsesRetainedEdge |= &updatedElement.edge( edge ) == retained;
    BOOST_CHECK( elementUsesRetainedEdge );
}

/** @brief Neighbor measures count overlapping high-valence fans once. */
BOOST_AUTO_TEST_CASE( pointNeighborMeasuresDeduplicateOverlappingFans )
{
    constexpr SurfaceMesh::size_type fanSize = 128;
    auto mesh = makeMesh();
    addPoint( *mesh, 0, 0, 0 );
    for ( SurfaceMesh::size_type pointId = 1; pointId <= fanSize; ++pointId )
    {
        double angle = 2 * std::acos( -1.0 ) * ( pointId - 1 ) / fanSize;
        addPoint( *mesh, pointId, std::cos( angle ), std::sin( angle ) );
    }
    auto const& first = mesh->point( 1 );
    auto const& second = mesh->point( 2 );
    addPoint( *mesh, fanSize + 1, 1.3 * ( first( 0 ) + second( 0 ) ) / 2,
              1.3 * ( first( 1 ) + second( 1 ) ) / 2 );

    auto addTriangle = [&mesh]( SurfaceMesh::size_type id, std::array<SurfaceMesh::size_type, 3> pointIds )
    {
        SurfaceMesh::element_type element;
        element.setId( id );
        element.setProcessId( 0 );
        element.setProcessIdInPartition( 0 );
        for ( uint16_type vertex = 0; vertex < 3; ++vertex )
            element.setPoint( vertex, mesh->point( pointIds[vertex] ) );
        mesh->addElement( std::move( element ) );
    };
    for ( SurfaceMesh::size_type id = 0; id < fanSize; ++id )
        addTriangle( id, {0, id + 1, ( id + 1 ) % fanSize + 1} );
    addTriangle( fanSize, {1, 2, fanSize + 1} );

    mesh->components().reset();
    mesh->components().set( MESH_UPDATE_FACES_MINIMAL | MESH_ADD_ELEMENTS_INFO );
    mesh->updateForUse();

    BOOST_CHECK_EQUAL( mesh->element( 0 ).pointElementNeighborIds().size(), fanSize + 1 );
    BOOST_CHECK_EQUAL( mesh->element( 2 ).pointElementNeighborIds().size(), fanSize );
    for ( auto const& [id, element] : mesh->elements() )
    {
        double expected = 0;
        for ( auto neighborId : element.pointElementNeighborIds() )
            expected += mesh->element( neighborId ).measure();
        BOOST_CHECK_CLOSE( element.measurePointElementNeighbors(), expected, 1e-8 );
    }
}
