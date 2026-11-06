/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file test_nedelec_interpolation.hpp
 * @brief Shared lowest-order Nedelec edge-interpolation regression.
 */
#pragma once

#include <feel/feelvf/vf.hpp>

#include <functional>
#include <string>
#include <vector>

namespace Feel::test
{
/**
 * @brief Return the canonical global orientation of a one-edge range.
 *
 * Feel++ orients global H(curl) edge DoFs from the smaller stable point
 * identifier to the larger one. Imported boundary entities may retain the
 * opposite Gmsh orientation, so their DSEL tangent must be corrected before
 * it is compared with a FunctionSpace coefficient.
 *
 * @tparam EdgeRangeType marked edge range type
 * @param range range containing one geometric edge
 * @return `+1` when the range follows global orientation, `-1` otherwise
 */
template <typename EdgeRangeType>
double canonicalEdgeOrientation( EdgeRangeType const& range )
{
    BOOST_REQUIRE_EQUAL( range.size(), 1 );
    auto const& entity = range.front();
    return entity.point( 0 ).id() < entity.point( 1 ).id()
               ? 1.0
               : -1.0;
}

/**
 * @brief Compare geometric edge moments with cell and boundary interpolation.
 *
 * The lowest-order NED1 dual consists of oriented tangential edge moments.
 * This helper evaluates those mathematical functionals directly with the
 * Geomap tangent, then compares their finite-element field with the public
 * `Element::on()` interpolation paths. The same kernel is used for boundary
 * edges of triangles and tetrahedra.
 *
 * @tparam MeshPtrType mesh smart-pointer type
 * @tparam SpacePtrType H(curl) function-space smart-pointer type
 * @tparam ExpressionType vector expression type
 * @tparam EdgeRangeFactory callable returning the marked geometric edge range
 * @tparam OrientationFactory callable returning the range-to-DoF sign
 * @param mesh one-cell simplex mesh
 * @param space lowest-order NED1 function space
 * @param edgeMarkers marker names in local Nedelec DoF order
 * @param expression vector field to interpolate
 * @param edgeRange marked-edge range factory
 * @param orientation range-to-global-DoF orientation policy
 * @param tolerance absolute coefficient and L2 tolerance
 */
template <typename MeshPtrType, typename SpacePtrType,
          typename ExpressionType, typename EdgeRangeFactory,
          typename OrientationFactory>
void checkLowestOrderNedelecInterpolation(
    MeshPtrType const& mesh,
    SpacePtrType const& space,
    std::vector<std::string> const& edgeMarkers,
    ExpressionType const& expression,
    EdgeRangeFactory&& edgeRange,
    OrientationFactory&& orientation,
    double tolerance = 1e-12 )
{
    BOOST_REQUIRE_EQUAL(
        static_cast<size_type>( space->nLocalDof() ),
        edgeMarkers.size() );

    auto edgeMoments = space->element();
    auto cellInterpolant = space->element();
    auto boundaryInterpolant = space->element();

    for ( size_type localDof = 0;
          localDof < edgeMarkers.size(); ++localDof )
    {
        auto const range =
            std::invoke( edgeRange, edgeMarkers[localDof] );
        auto const orientationSign =
            std::invoke( orientation, range );
        auto const moment = integrate(
            _range = range,
            _expr = cst( orientationSign ) *
                    trans( T() ) * expression )
                                .evaluate()( 0, 0 );
        edgeMoments( localDof ) = moment;
    }

    cellInterpolant.on(
        _range = elements( mesh ), _expr = expression );
    boundaryInterpolant.on(
        _range = boundaryfaces( mesh ), _expr = expression );

    for ( size_type localDof = 0;
          localDof < edgeMarkers.size(); ++localDof )
    {
        BOOST_CHECK_SMALL(
            edgeMoments( localDof ) -
                cellInterpolant( localDof ),
            tolerance );
        BOOST_CHECK_SMALL(
            boundaryInterpolant( localDof ) -
                cellInterpolant( localDof ),
            tolerance );
    }

    auto const cellError = normL2(
        _range = elements( mesh ),
        _expr = idv( edgeMoments ) -
                idv( cellInterpolant ) );
    auto const boundaryError = normL2(
        _range = elements( mesh ),
        _expr = idv( boundaryInterpolant ) -
                idv( cellInterpolant ) );

    BOOST_TEST_MESSAGE(
        "NED1 edge-moment/cell L2 error = " << cellError );
    BOOST_TEST_MESSAGE(
        "NED1 boundary/cell L2 error = " << boundaryError );
    BOOST_CHECK_SMALL( cellError, tolerance );
    BOOST_CHECK_SMALL( boundaryError, tolerance );
}
} // namespace Feel::test
