/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_geomap_frames.cpp
 * @brief Public Geomap facet-frame and tangent-projector tests.
 */
#define BOOST_TEST_MODULE test_geomap_frames
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitcircle.hpp>
#include <feel/feelfilters/unithypercube.hpp>
#include <feel/feelpoly/geomap.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>
#include <set>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/** Detect contexts for which a facet has one mathematically unique tangent. */
template<typename ContextType>
concept single_tangent_context = requires( ContextType const& context )
{
    context.unitTangent( 0 );
};

using mesh2_type = Mesh<Simplex<2,2>>;
using mesh3_type = Mesh<Simplex<3,1>>;
using context2_type = typename mesh2_type::gm_type::template Context<
    typename mesh2_type::element_type,1>;
using context3_type = typename mesh3_type::gm_type::template Context<
    typename mesh3_type::element_type,1>;

static_assert( single_tangent_context<context2_type> );
static_assert( !single_tangent_context<context3_type> );
} // namespace

BOOST_AUTO_TEST_SUITE( geomap_frames_suite )

/**
 * Verify the public frame contracts on every face of curved quadratic cells.
 */
BOOST_AUTO_TEST_CASE( curved_facet_frames_and_projectors )
{
    if ( Environment::worldComm().localSize() > 1 )
    {
        BOOST_TEST_MESSAGE( "curved Geomap frame contract is serial-only" );
        return;
    }

    auto mesh = unitCircle<2>( 0.35 );
    BOOST_REQUIRE( mesh );
    auto gm = mesh->gm();
    IMGeneral<2,double,Simplex> quadrature( 4 );
    auto facePrecompute = gm->preComputeOnFaces(
        gm, quadrature.allfpoints() );
    constexpr size_type context = vm::POINT | vm::JACOBIAN | vm::KB |
                                  vm::NORMAL | vm::TANGENT |
                                  vm::LOCAL_BASIS;

    std::set<uint16_type> permutations;
    double maximumFrameVariation = 0.0;
    bool initialized = false;
    Eigen::Matrix2d firstPiolaFrame = Eigen::Matrix2d::Zero();

    for ( auto const& elementWrap : elements( mesh ) )
    {
        auto const& element = unwrap_ref( elementWrap );
        for ( uint16_type face = 0;
              face < mesh2_type::element_type::numTopologicalFaces; ++face )
        {
            auto gmc = gm->template context<context>(
                element, facePrecompute, face );
            permutations.insert( gmc->permutation().value() );
            for ( int q = 0; q < gmc->nPoints(); ++q )
            {
                auto const& frame = gmc->orthonormalFrame( q );
                auto const tangents = gmc->orientedTangentFrame( q );
                auto const& projector = gmc->tangentProjector( q );
                auto const piolaFrame = gmc->piolaTraceFrame( face, q );

                BOOST_CHECK_SMALL(
                    ( frame.transpose()*frame-
                      Eigen::Matrix2d::Identity() ).norm(), 3e-12 );
                BOOST_CHECK_GT( frame.determinant(), 0.0 );
                BOOST_CHECK_GT(
                    frame.col( 0 ).dot( gmc->unitNormal( q ) ), 1.0-3e-12 );
                BOOST_CHECK_SMALL(
                    ( projector-projector.transpose() ).norm(), 3e-12 );
                BOOST_CHECK_SMALL(
                    ( projector*projector-projector ).norm(), 3e-12 );
                BOOST_CHECK_SMALL(
                    ( projector*gmc->unitNormal( q ) ).norm(), 3e-12 );
                BOOST_CHECK_SMALL(
                    ( projector-geometry::subspaceProjector( tangents ) ).norm(),
                    3e-12 );
                BOOST_CHECK_SMALL(
                    ( piolaFrame.col( 0 )-
                      gmc->piolaContravariantFacetNormal( face, q ) ).norm(),
                    3e-12 );

                auto const edgeTangent =
                    gmc->piolaCovariantEdgeTangent( face, q );
                BOOST_CHECK_SMALL(
                    std::abs( piolaFrame.col( 1 ).dot( edgeTangent ) )-
                        piolaFrame.col( 1 ).norm()*edgeTangent.norm(),
                    3e-12 );

                if ( !initialized )
                {
                    firstPiolaFrame = piolaFrame;
                    initialized = true;
                }
                else
                    maximumFrameVariation = std::max(
                        maximumFrameVariation,
                        ( piolaFrame-firstPiolaFrame ).norm() );
            }
        }
    }

    BOOST_CHECK( initialized );
    BOOST_CHECK_GE( permutations.size(), 2 );
    BOOST_CHECK_GT( maximumFrameVariation, 1e-5 );
}

/**
 * Verify three-dimensional simplex and tensor-product facet-frame contracts.
 */
BOOST_AUTO_TEST_CASE( three_dimensional_facet_frames )
{
    if ( Environment::worldComm().localSize() > 1 )
    {
        BOOST_TEST_MESSAGE( "3D Geomap frame contract is serial-only" );
        return;
    }

    auto checkFrames = []<typename MeshType, typename QuadratureType>(
                           std::shared_ptr<MeshType> const& mesh,
                           QuadratureType const& quadrature,
                           bool tensorProduct,
                           bool requireInteriorInterface )
    {
        using element_type = typename MeshType::element_type;
        BOOST_REQUIRE( mesh );
        auto gm = mesh->gm();
        auto facePrecompute = gm->preComputeOnFaces(
            gm, quadrature.allfpoints() );
        constexpr size_type context = vm::POINT | vm::JACOBIAN | vm::KB |
                                      vm::NORMAL | vm::TANGENT |
                                      vm::LOCAL_BASIS;
        std::set<uint16_type> permutations;
        std::size_t checkedFaces = 0;

        for ( auto const& elementWrap : elements( mesh ) )
        {
            auto const& element = unwrap_ref( elementWrap );
            for ( uint16_type face = 0;
                  face < element_type::numTopologicalFaces; ++face )
            {
                auto gmc = gm->template context<context>(
                    element, facePrecompute, face );
                permutations.insert( gmc->permutation().value() );
                auto const referenceTangents =
                    gmc->referenceFacetTangentFrame( face );
                Eigen::Matrix3d referenceFrame;
                referenceFrame.col( 0 ) =
                    gm->referenceConvex().normal( face );
                referenceFrame.rightCols( 2 ) = referenceTangents;
                BOOST_CHECK_GT( referenceFrame.determinant(), 0.0 );
                BOOST_CHECK_GT( gmc->referenceFacetMeasure( face ), 0.0 );
                if ( tensorProduct )
                    BOOST_CHECK_CLOSE(
                        gmc->referenceFacetMeasure( face ), 4.0, 1e-12 );

                for ( int q = 0; q < gmc->nPoints(); ++q )
                {
                    auto const& frame = gmc->orthonormalFrame( q );
                    auto const tangents = gmc->orientedTangentFrame( q );
                    auto const& projector = gmc->tangentProjector( q );
                    auto const piolaFrame = gmc->piolaTraceFrame( face, q );

                    BOOST_CHECK_SMALL(
                        ( frame.transpose()*frame-
                          Eigen::Matrix3d::Identity() ).norm(), 5e-12 );
                    BOOST_CHECK_GT( frame.determinant(), 0.0 );
                    BOOST_CHECK_GT(
                        frame.col( 0 ).dot( gmc->unitNormal( q ) ),
                        1.0-5e-12 );
                    BOOST_CHECK_SMALL(
                        ( projector-projector.transpose() ).norm(), 5e-12 );
                    BOOST_CHECK_SMALL(
                        ( projector*projector-projector ).norm(), 5e-12 );
                    BOOST_CHECK_SMALL(
                        ( projector*gmc->unitNormal( q ) ).norm(), 5e-12 );
                    BOOST_CHECK_SMALL(
                        ( projector-
                          geometry::subspaceProjector( tangents ) ).norm(),
                        5e-12 );
                    BOOST_CHECK_SMALL(
                        ( piolaFrame.col( 0 )-
                          gmc->piolaContravariantFacetNormal( face, q ) ).norm(),
                        5e-12 );
                    BOOST_CHECK_SMALL(
                        ( piolaFrame.rightCols( 2 )-
                          gmc->K( q )*referenceTangents ).norm(), 5e-12 );
                }
                ++checkedFaces;
            }
        }

        BOOST_CHECK_GT( checkedFaces, 0 );
        BOOST_CHECK( !permutations.empty() );

        std::size_t checkedInterfaces = 0;
        for ( auto faceIt = mesh->beginFace(); faceIt != mesh->endFace();
              ++faceIt )
        {
            auto const& interface = faceIt->second;
            if ( interface.isOnBoundary() || !interface.isConnectedTo0() ||
                 !interface.isConnectedTo1() )
                continue;

            auto gmc0 = gm->template context<context>(
                mesh->element( interface.ad_first() ), facePrecompute,
                interface.pos_first() );
            auto gmc1 = gm->template context<context>(
                mesh->element( interface.ad_second() ), facePrecompute,
                interface.pos_second() );

            for ( int q0 = 0; q0 < gmc0->nPoints(); ++q0 )
            {
                int matchingPoint = -1;
                double minimumSquaredDistance =
                    std::numeric_limits<double>::infinity();
                for ( int q1 = 0; q1 < gmc1->nPoints(); ++q1 )
                {
                    double squaredDistance = 0.0;
                    for ( int component = 0; component < 3; ++component )
                    {
                        double const difference =
                            gmc0->xReal( q0 )( component )-
                            gmc1->xReal( q1 )( component );
                        squaredDistance += difference*difference;
                    }
                    if ( squaredDistance < minimumSquaredDistance )
                    {
                        minimumSquaredDistance = squaredDistance;
                        matchingPoint = q1;
                    }
                }

                BOOST_REQUIRE_GE( matchingPoint, 0 );
                BOOST_CHECK_SMALL( minimumSquaredDistance, 2e-24 );
                BOOST_CHECK_SMALL(
                    gmc0->unitNormal( q0 ).dot(
                        gmc1->unitNormal( matchingPoint ) )+1.0,
                    5e-12 );
                BOOST_CHECK_SMALL(
                    ( gmc0->tangentProjector( q0 )-
                      gmc1->tangentProjector( matchingPoint ) ).norm(),
                    5e-12 );
            }
            ++checkedInterfaces;
        }
        if ( requireInteriorInterface )
            BOOST_CHECK_GT( checkedInterfaces, 0 );
    };

    auto simplexMesh = unitCube( 0.8 );
    IMGeneral<3,double,Simplex> simplexQuadrature( 3 );
    checkFrames( simplexMesh, simplexQuadrature, false, true );

    auto tensorProductMesh = unitHypercube<3,Hypercube<3>>( 0.8 );
    IMGeneral<3,double,Hypercube> tensorProductQuadrature( 3 );
    checkFrames( tensorProductMesh, tensorProductQuadrature, true, false );
}

BOOST_AUTO_TEST_SUITE_END()
