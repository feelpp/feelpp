/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_localization.cpp
 * @brief C++20 Localization API, policy, and compatibility tests.
 */
#define BOOST_TEST_MODULE localization testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/unitsquare.hpp>

#include <type_traits>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
/**
 * @brief Check named and compatibility localization APIs on one mesh.
 * @tparam MeshType mesh type
 * @param mesh distributed simplex mesh
 */
template <typename MeshType>
void checkLocalizationApi( std::shared_ptr<MeshType> const& mesh )
{
    using localization_type = Localization<MeshType>;
    using size_type = typename localization_type::size_type;
    using matrix_node_type = typename localization_type::matrix_node_type;

    auto elementRange = elements( mesh );
    auto elementIt = elementRange.begin();
    if ( elementIt == elementRange.end() )
    {
        BOOST_TEST_MESSAGE( "rank owns no element; localization check skipped" );
        return;
    }

    auto const& element = unwrap_ref( *elementIt );
    auto const physicalPoint = element.barycenter();
    auto localization = mesh->tool_localization();
    BOOST_REQUIRE( localization );
    localization->setExtrapolation( false );
    localization->updateForUse();

    auto pointLocation =
        localization->locateInElement( element.id(), physicalPoint );
    BOOST_REQUIRE( pointLocation.found );
    BOOST_CHECK_EQUAL( pointLocation.referencePoint.size(), MeshType::nDim );

    auto elementLocation = localization->locateElement( physicalPoint );
    BOOST_REQUIRE( elementLocation.found );
    BOOST_CHECK_EQUAL( elementLocation.elementId, element.id() );

    matrix_node_type points( MeshType::nRealDim, 1 );
    ublas::column( points, 0 ) = physicalPoint;
    auto unconstrained = localization->localize(
        points, invalid_v<size_type> );
    BOOST_REQUIRE_EQUAL( unconstrained.found.size(), 1 );
    BOOST_CHECK( unconstrained.found[0] );
    BOOST_CHECK_EQUAL( unconstrained.lastElementId, element.id() );
    BOOST_REQUIRE_EQUAL( localization->resultAnalysis().size(), 1 );
    BOOST_CHECK_EQUAL(
        localization->resultAnalysis().begin()->first, element.id() );

    auto conforming = localization->template localize<
        LocalizationPolicy::RequireContainedPoints>(
            points, invalid_v<size_type>, element.vertices() );
    BOOST_REQUIRE_EQUAL( conforming.found.size(), 1 );
    BOOST_CHECK( conforming.found[0] );
    BOOST_CHECK_EQUAL( conforming.lastElementId, element.id() );

    auto conformingBoolean = localization->localize(
        points, invalid_v<size_type>, element.vertices(), std::true_type{} );
    BOOST_REQUIRE_EQUAL( conformingBoolean.found.size(), 1 );
    BOOST_CHECK( conformingBoolean.found[0] );
    BOOST_CHECK_EQUAL(
        conformingBoolean.lastElementId, conforming.lastElementId );

    matrix_node_type invalidRequiredPoints = element.vertices();
    for ( size_type component = 0;
          component < invalidRequiredPoints.size1();
          ++component )
        invalidRequiredPoints( component, 0 ) += 10.0;
    auto rejected = localization->template locateInElement<
        LocalizationPolicy::RequireContainedPoints>(
            element.id(), physicalPoint, invalidRequiredPoints );
    BOOST_CHECK( !rejected.found );

    auto legacy = localization->isIn( element.id(), physicalPoint );
    BOOST_CHECK_EQUAL( boost::get<0>( legacy ), pointLocation.found );
    BOOST_CHECK_SMALL(
        ublas::norm_2(
            boost::get<1>( legacy ) - pointLocation.referencePoint ),
        1e-12 );
    BOOST_CHECK_SMALL(
        std::abs( boost::get<2>( legacy ) - pointLocation.distance ),
        1e-12 );

    auto outsidePoint = physicalPoint;
    for ( size_type component = 0; component < outsidePoint.size(); ++component )
        outsidePoint( component ) += 10.0;
    auto outside = localization->locateElement( outsidePoint );
    BOOST_CHECK( !outside.found );
}
} // namespace

BOOST_AUTO_TEST_SUITE( localization_suite )

/** @test Exercise localization policies and named results on triangles. */
BOOST_AUTO_TEST_CASE( simplex_2d )
{
    checkLocalizationApi( unitSquare( 0.6 ) );
}

/** @test Exercise localization policies and named results on tetrahedra. */
BOOST_AUTO_TEST_CASE( simplex_3d )
{
    checkLocalizationApi( unitCube( 0.9 ) );
}

BOOST_AUTO_TEST_SUITE_END()
