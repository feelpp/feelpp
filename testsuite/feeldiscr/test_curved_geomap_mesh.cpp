/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_curved_geomap_mesh.cpp
 * @brief Imported curved-mesh validity and shared-node contracts.
 */

#define BOOST_TEST_MODULE curved geomap mesh testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelfilters/unitsphere.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelpoly/im.hpp>

#include <cmath>
#include <set>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/**
 * @brief Check pointwise map validity in cell interiors and on cell facets.
 *
 * @tparam Dim topological and ambient dimension
 * @tparam GeometryOrder geometric polynomial order
 * @param mesh curved simplex mesh
 */
template <int Dim, int GeometryOrder>
void checkMappingValidity(
    std::shared_ptr<Mesh<Simplex<Dim, GeometryOrder>>> const& mesh )
{
    using mesh_type = Mesh<Simplex<Dim, GeometryOrder>>;
    using element_type = typename mesh_type::element_type;
    constexpr size_type context =
        vm::POINT | vm::JACOBIAN | vm::KB | vm::SECOND_DERIVATIVE;

    auto gm = mesh->gm();
    IMGeneral<Dim, double, Simplex> quadrature( 14 );
    auto precompute = gm->preCompute( gm, quadrature.points() );
    auto facePrecompute =
        gm->preComputeOnFaces( gm, quadrature.allfpoints() );

    for ( auto const& elementWrap : elements( mesh ) )
    {
        auto const& element = unwrap_ref( elementWrap );
        auto gmc = gm->template context<context>( element, precompute );
        for ( int q = 0; q < gmc->nPoints(); ++q )
        {
            BOOST_CHECK( gmc->isValid( q ) );
            BOOST_CHECK( std::isfinite( gmc->conditionNumber( q ) ) );
            BOOST_CHECK_GT( gmc->signedJacobian( q ), 0.0 );
        }

        for ( uint16_type face = 0;
              face < element_type::numTopologicalFaces; ++face )
        {
            auto faceContext = gm->template context<context>(
                element, facePrecompute, face );
            for ( int q = 0; q < faceContext->nPoints(); ++q )
            {
                BOOST_CHECK( faceContext->isValid( q ) );
                BOOST_CHECK_GT(
                    faceContext->signedJacobian( q ), 0.0 );
            }
        }
    }
}

/**
 * @brief Check canonical face-interior node ownership on tetrahedra.
 *
 * @tparam GeometryOrder geometric polynomial order
 * @param mesh curved tetrahedral mesh
 */
template <int GeometryOrder>
void checkSharedFaceInteriorNodes(
    std::shared_ptr<Mesh<Simplex<3, GeometryOrder>>> const& mesh )
{
    if constexpr ( GeometryOrder < 3 )
        return;
    else
    {
        using mesh_type = Mesh<Simplex<3, GeometryOrder>>;
        using element_type = typename mesh_type::element_type;
        using face_type = typename mesh_type::face_type;
        using size_type = typename mesh_type::size_type;
        constexpr uint16_type firstFaceInteriorPoint =
            face_type::numVertices * face_type::nbPtsPerVertex +
            face_type::numEdges * face_type::nbPtsPerEdge;

        std::size_t checkedFaces = 0;
        for ( auto const& faceWrap : internalfaces( mesh ) )
        {
            auto const& face = unwrap_ref( faceWrap );
            BOOST_REQUIRE( face.isConnectedTo0() );
            BOOST_REQUIRE( face.isConnectedTo1() );

            std::set<size_type> canonicalIds;
            std::set<size_type> firstElementIds;
            std::set<size_type> secondElementIds;
            for ( uint16_type point = firstFaceInteriorPoint;
                  point < face_type::numPoints; ++point )
            {
                canonicalIds.insert( face.point( point ).id() );
                firstElementIds.insert(
                    face.element0()
                        .point( element_type::fToP(
                            face.pos_first(), point ) )
                        .id() );
                secondElementIds.insert(
                    face.element1()
                        .point( element_type::fToP(
                            face.pos_second(), point ) )
                        .id() );
            }

            BOOST_CHECK_EQUAL_COLLECTIONS(
                canonicalIds.begin(), canonicalIds.end(),
                firstElementIds.begin(), firstElementIds.end() );
            BOOST_CHECK_EQUAL_COLLECTIONS(
                canonicalIds.begin(), canonicalIds.end(),
                secondElementIds.begin(), secondElementIds.end() );
            ++checkedFaces;
        }
        BOOST_CHECK_GT( checkedFaces, 0 );
    }
}

/**
 * @brief Exercise one curved geometric order in two and three dimensions.
 *
 * @tparam GeometryOrder geometric polynomial order
 */
template <int GeometryOrder>
void checkGeometryOrder()
{
    auto circle = unitCircle<GeometryOrder>( 0.55 );
    auto sphere = unitSphere<GeometryOrder>( 0.9 );
    BOOST_REQUIRE( circle );
    BOOST_REQUIRE( sphere );

    checkMappingValidity<2, GeometryOrder>( circle );
    checkMappingValidity<3, GeometryOrder>( sphere );
    checkSharedFaceInteriorNodes<GeometryOrder>( sphere );
}
} // namespace

BOOST_AUTO_TEST_SUITE( curved_geomap_mesh )

/** Validate imported quadratic circle and sphere mappings. */
BOOST_AUTO_TEST_CASE( quadratic_geometry )
{
    checkGeometryOrder<2>();
}

/** Validate imported cubic circle and sphere mappings and shared face nodes. */
BOOST_AUTO_TEST_CASE( cubic_geometry )
{
#if FEELPP_MESH_MAX_ORDER >= 3
    checkGeometryOrder<3>();
#else
    BOOST_TEST_MESSAGE(
        "cubic mesh validation requires FEELPP_MESH_MAX_ORDER >= 3" );
#endif
}

BOOST_AUTO_TEST_SUITE_END()
