/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_geoentity_views.cpp
 * @brief Validate mathematical mesh-entity views over compact storage.
 */

#define BOOST_TEST_MAIN
#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelmesh/concepts.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <ranges>
#include <sstream>
#include <type_traits>
#include <vector>

FEELPP_ENVIRONMENT_NO_OPTIONS

using namespace Feel;

namespace
{
using test_mesh_type = Mesh<Simplex<2, 2>>;

static_assert( GeometricPointViewType<
               test_mesh_type::point_view_type> );
static_assert( GeometricEntityViewType<
               test_mesh_type::edge_view_type> );
static_assert( GeometricEdgeViewType<
               test_mesh_type::edge_view_type> );
static_assert( GeometricEntityViewType<
               test_mesh_type::face_view_type> );
static_assert( GeometricFaceViewType<
               test_mesh_type::face_view_type> );
static_assert( GeometricEntityViewType<
               test_mesh_type::element_view_type> );
static_assert( GeometricCellViewType<
               test_mesh_type::element_view_type> );
static_assert( std::is_trivially_copyable_v<
               test_mesh_type::element_view_type> );
} // namespace

/**
 * @brief Check canonical connectivity and zero-copy point-coordinate views.
 */
BOOST_AUTO_TEST_CASE( entity_view_matches_compatibility_geometry )
{
    auto mesh = createGMSHMesh(
        _mesh = new test_mesh_type,
        _desc = domain(
            _name = "entity-view-geometry",
            _shape = "hypercube",
            _dim = 2, _h = 0.35 ),
        _verbose = 0 );

    auto const& element = mesh->beginElement()->second;
    auto const view = mesh->elementView( element.id() );
    auto const copiedView = view;

    BOOST_CHECK( view == copiedView );
    BOOST_CHECK_EQUAL( view.id(), element.id() );
    BOOST_CHECK_EQUAL( view.pointCount(), element.nPoints() );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        view.pointIds().begin(), view.pointIds().end(),
        element.pointIds().begin(), element.pointIds().end() );

    auto const localIndices = view.pointLocalIndices();
    auto const allCoordinates =
        static_cast<test_mesh_type const&>( *mesh )
            .pointCoordinates();
    auto const viewCoordinates = view.coordinates();
    auto const legacyCoordinates = element.G();
    BOOST_REQUIRE_EQUAL(
        viewCoordinates.cols(), legacyCoordinates.size2() );
    for ( uint16_type localPoint = 0;
          localPoint < view.pointCount(); ++localPoint )
    {
        auto const point = view.point( localPoint );
        BOOST_CHECK_EQUAL(
            point.id(), view.pointIds()[localPoint] );
        BOOST_CHECK_EQUAL(
            point.localIndex(), localIndices[localPoint] );
        for ( uint16_type component = 0;
              component < test_mesh_type::nRealDim; ++component )
        {
            BOOST_CHECK_CLOSE(
                point.coordinate( component ),
                allCoordinates(
                    component, localIndices[localPoint] ),
                1e-12 );
            BOOST_CHECK_CLOSE(
                viewCoordinates( component, localPoint ),
                legacyCoordinates( component, localPoint ),
                1e-12 );
        }
    }
}

/**
 * @brief Check that stored entities do not retain point-pointer ownership.
 *
 * Compact connectivity is canonical after mesh insertion. Direct `point()`
 * access resolves through the mesh-wide owner without allocating a pointer
 * vector; the legacy `points()` facade materializes one only on demand.
 */
BOOST_AUTO_TEST_CASE( stored_entities_use_lazy_point_facades )
{
    auto mesh = createGMSHMesh(
        _mesh = new test_mesh_type,
        _desc = domain(
            _name = "entity-view-lazy-point-facade",
            _shape = "hypercube",
            _dim = 2, _h = 0.3 ),
        _verbose = 0 );

    for ( auto const& [id, element] : mesh->elements() )
    {
        BOOST_TEST_CONTEXT( "element " << id )
        {
            BOOST_CHECK( !element.hasPointPointerCache() );
            for ( uint16_type localPoint = 0;
                  localPoint < element.nPointsRuntime();
                  ++localPoint )
                BOOST_CHECK_EQUAL(
                    element.point( localPoint ).id(),
                    element.pointIds()[localPoint] );
            BOOST_CHECK( !element.hasPointPointerCache() );
        }
    }
    for ( auto const& [id, face] : mesh->faces() )
    {
        BOOST_TEST_CONTEXT( "face " << id )
        {
            BOOST_CHECK( !face.hasPointPointerCache() );
            BOOST_CHECK_EQUAL(
                face.point( 0 ).id(), face.pointIds()[0] );
            BOOST_CHECK( !face.hasPointPointerCache() );
        }
    }

    auto const& firstElement =
        mesh->beginElement()->second;
    auto const& compatibilityPoints =
        firstElement.points();
    BOOST_CHECK( firstElement.hasPointPointerCache() );
    BOOST_CHECK_EQUAL(
        compatibilityPoints.size(),
        firstElement.pointIds().size() );
}

/**
 * @brief Check cell subentity, adjacency, and lazy range views.
 */
BOOST_AUTO_TEST_CASE( entity_view_traverses_mathematical_topology )
{
    auto mesh = createGMSHMesh(
        _mesh = new test_mesh_type,
        _desc = domain(
            _name = "entity-view-topology",
            _shape = "hypercube",
            _dim = 2, _h = 0.25 ),
        _verbose = 0 );

    size_type viewedElements = 0;
    bool foundNeighbor = false;
    for ( auto const cell : elementViews( mesh ) )
    {
        ++viewedElements;
        BOOST_CHECK_EQUAL(
            cell.entity().id(), cell.id() );

        for ( uint16_type localFace = 0;
              localFace < cell.entity().numTopologicalFaces;
              ++localFace )
        {
            auto const face = cell.face( localFace );
            BOOST_CHECK_EQUAL(
                face.id(),
                cell.entity().face( localFace ).id() );
            BOOST_CHECK_EQUAL(
                face.pointCount(),
                cell.entity().face( localFace ).nPoints() );

            auto const neighbor = cell.neighbor( localFace );
            if ( neighbor )
            {
                foundNeighbor = true;
                BOOST_CHECK_EQUAL(
                    neighbor->id(),
                    cell.entity().neighbor( localFace ) );
            }
        }
    }

    BOOST_CHECK_EQUAL( viewedElements, mesh->numElements() );
    BOOST_CHECK( foundNeighbor );
    BOOST_CHECK_EQUAL(
        std::ranges::distance( pointViews( mesh ) ),
        mesh->numPoints() );
    BOOST_CHECK_EQUAL(
        std::ranges::distance( faceViews( mesh ) ),
        mesh->numFaces() );
}

/**
 * @brief Check Context and Inverse parity through compact entity views.
 */
BOOST_AUTO_TEST_CASE( entity_view_drives_geometric_mapping )
{
    auto mesh = createGMSHMesh(
        _mesh = new test_mesh_type,
        _desc = domain(
            _name = "entity-view-geomap",
            _shape = "hypercube",
            _dim = 2, _h = 0.3 ),
        _verbose = 0 );

    auto elementIterator = mesh->beginElement();
    ++elementIterator;
    BOOST_REQUIRE( elementIterator != mesh->endElement() );
    auto const& secondElement = elementIterator->second;
    auto const secondView =
        mesh->elementView( secondElement.id() );

    auto const geometricMapping = mesh->gm();
    typename test_mesh_type::reference_convex_type
        referenceElement;
    auto const precompute =
        geometricMapping->preCompute( referenceElement.points() );
    constexpr size_type context =
        vm::POINT | vm::JACOBIAN | vm::KB;

    auto legacyContext =
        geometricMapping->template context<context>(
            secondElement, precompute );
    auto viewContext =
        geometricMapping->template context<context>(
            secondView, precompute );

    BOOST_REQUIRE_EQUAL(
        legacyContext->nPoints(), viewContext->nPoints() );
    for ( uint16_type q = 0;
          q < legacyContext->nPoints(); ++q )
    {
        BOOST_CHECK_SMALL(
            ublas::norm_2(
                legacyContext->xReal( q ) -
                viewContext->xReal( q ) ),
            1e-12 );
        BOOST_CHECK_CLOSE_FRACTION(
            legacyContext->J( q ),
            viewContext->J( q ), 1e-12 );
    }

    using inverse_type =
        typename test_mesh_type::gm_type::Inverse;
    inverse_type legacyInverse(
        geometricMapping, secondElement );
    inverse_type viewInverse(
        geometricMapping, secondView );
    auto const physicalPoint =
        legacyContext->xReal( 0 );
    legacyInverse.setXReal( physicalPoint );
    viewInverse.setXReal( physicalPoint );
    BOOST_CHECK_SMALL(
        ublas::norm_2(
            legacyInverse.xRef() - viewInverse.xRef() ),
        1e-12 );
}

/**
 * @brief Check that DoF lookup accepts the common entity-view contract.
 */
BOOST_AUTO_TEST_CASE( entity_view_drives_dof_lookup )
{
    auto mesh = createGMSHMesh(
        _mesh = new test_mesh_type,
        _desc = domain(
            _name = "entity-view-dof",
            _shape = "hypercube",
            _dim = 2, _h = 0.3 ),
        _verbose = 0 );
    auto const space = Pch<1>( mesh );
    auto const dof = space->dof();
    auto const& element = mesh->beginElement()->second;
    auto const view = mesh->elementView( element.id() );

    auto const localDofCount =
        static_cast<uint16_type>(
            dof->nLocalDof( true ) );
    for ( uint16_type localDof = 0;
          localDof < localDofCount; ++localDof )
    {
        auto const fromElement =
            dof->localToGlobal( element, localDof );
        auto const fromView =
            dof->localToGlobal( view, localDof );
        BOOST_CHECK_EQUAL(
            fromElement.index(), fromView.index() );
        BOOST_CHECK_EQUAL(
            dof->localDofInElement( element, localDof ),
            dof->localDofInElement( view, localDof ) );
    }
}

/**
 * @brief Check that views are reconstructed after mesh serialization.
 *
 * Views deliberately serialize no owner pointer. A restored mesh creates new
 * value-semantic handles over its decoded compact connectivity and coordinate
 * storage.
 */
BOOST_AUTO_TEST_CASE( entity_views_rebind_after_mesh_restart )
{
    if ( Environment::worldComm().globalSize() != 1 )
    {
        BOOST_TEST_MESSAGE(
            "full-mesh archive restart is validated in serial; "
            "distributed view behavior is covered by the other cases" );
        return;
    }

    auto mesh = createGMSHMesh(
        _mesh = new test_mesh_type,
        _desc = domain(
            _name = "entity-view-restart",
            _shape = "hypercube",
            _dim = 2, _h = 0.35 ),
        _verbose = 0 );
    auto const sourceElementId =
        mesh->beginElement()->second.id();
    auto const sourcePointIds =
        mesh->elementView( sourceElementId ).pointIds();
    std::vector<test_mesh_type::index_type> expectedPointIds(
        sourcePointIds.begin(), sourcePointIds.end() );

    std::stringstream stream(
        std::ios::in | std::ios::out | std::ios::binary );
    {
        boost::archive::binary_oarchive archive( stream );
        archive << *mesh;
    }

    test_mesh_type restored;
    {
        boost::archive::binary_iarchive archive( stream );
        archive >> restored;
    }

    auto const restoredView =
        restored.elementView( sourceElementId );
    auto const& restoredElement =
        restored.element( sourceElementId );
    BOOST_CHECK( !restoredElement.hasPointPointerCache() );
    BOOST_CHECK_EQUAL_COLLECTIONS(
        restoredView.pointIds().begin(),
        restoredView.pointIds().end(),
        expectedPointIds.begin(), expectedPointIds.end() );
    auto const restoredCoordinates =
        restoredView.coordinates();
    BOOST_CHECK_EQUAL(
        restoredCoordinates.cols(),
        restoredView.pointCount() );
    BOOST_CHECK_EQUAL(
        restoredElement.nNeighbors(),
        restoredElement.numNeighbors );
    for ( uint16_type localFace = 0;
          localFace < restoredElement.nNeighbors();
          ++localFace )
    {
        auto const neighbor = restoredView.neighbor( localFace );
        BOOST_CHECK_EQUAL(
            static_cast<bool>( neighbor ),
            restoredElement.hasNeighbor( localFace ) );
    }
}
