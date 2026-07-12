/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file hcurlfunctionals.hpp @brief H(curl) family dual-functional builders. */

#ifndef FEELPP_FEELPOLY_HCURLFUNCTIONALS_HPP
#define FEELPP_FEELPOLY_HCURLFUNCTIONALS_HPP 1

#include <feel/feelpoly/tracefunctionals.hpp>

#include <algorithm>
#include <iterator>
#include <vector>

namespace Feel::detail
{

/**
 * Append the canonical edge-tangential point functionals of a reference
 * H(curl) element and populate its element/facet representative point tables.
 *
 * This is the common lowest-order/point-variant kernel. Higher-order moment
 * variants use the same entity traversal but supply polynomial weights in the
 * H(curl) moment layer.
 *
 * @tparam Space H(curl) primal polynomial space
 * @tparam ReferenceConvexType reference cell type
 * @param space primal space on which the edge functionals act
 * @param referenceConvex reference cell providing edge points and tangents
 * @param dofPerEdge number of ordered functionals attached to each edge
 * @param elementPoints destination table of representative element points
 * @param facetPoints destination tables for facet-closure points
 * @param functionals ordered destination dual block
 */
template<typename Space, typename ReferenceConvexType>
void
appendHCurlEdgeTangentialPointFunctionals(
    Space const& space,
    ReferenceConvexType const& referenceConvex,
    uint16_type dofPerEdge,
    typename Space::points_type& elementPoints,
    std::vector<typename Space::points_type>& facetPoints,
    std::vector<Functional<Space>>& functionals )
{
    using value_type = typename Space::value_type;
    using points_type = typename Space::points_type;
    using node_type = typename ReferenceConvexType::node_type;
    using convex_type = typename ReferenceConvexType::super;
    using face_type = typename convex_type::topological_face_type;

    static constexpr uint16_type nDim = Space::nDim;
    static_assert( nDim == 2 || nDim == 3,
                   "H(curl) simplex edge functionals currently support triangles and tetrahedra." );

    const uint16_type edgesPerFacet = nDim == 2
                                          ? face_type::numEdges
                                          : face_type::numTopologicalFaces;
    uint16_type pointOffset = 0;
    for ( int edge = referenceConvex.entityRange( 1 ).begin();
          edge < referenceConvex.entityRange( 1 ).end();
          ++edge )
    {
        points_type edgePoints( referenceConvex.makePoints( 1, edge ) );

        for ( int facet = referenceConvex.entityRange( nDim-1 ).begin();
              facet < referenceConvex.entityRange( nDim-1 ).end();
              ++facet )
        {
            for ( uint16_type localEdge = 0; localEdge < edgesPerFacet; ++localEdge )
            {
                const int facetEdge = nDim == 2
                                          ? convex_type::f2e( facet, facet )
                                          : convex_type::f2e( facet, localEdge );
                if ( facetEdge == edge )
                {
                    ublas::subrange( facetPoints[facet], 0, nDim,
                                     dofPerEdge*localEdge,
                                     dofPerEdge*( localEdge+1 ) ) = edgePoints;
                }
            }
        }

        if ( edgePoints.size2() != 0 )
        {
            ublas::subrange( elementPoints, 0, nDim,
                             pointOffset, pointOffset+edgePoints.size2() ) = edgePoints;
            pointOffset += static_cast<uint16_type>( edgePoints.size2() );
        }

        node_type tangent( nDim );
        em_node_type<value_type> mappedTangent( tangent.data().begin(), tangent.size() );
        mappedTangent = referenceConvex.tangent( edge );
        functional::appendDirectionalPointFunctionals( space, tangent, edgePoints, functionals );
    }
}

/**
 * Append the legacy weighted two-point edge functional used by the currently
 * supported Nedelec-second-kind triangle. It is isolated here so the NED2
 * family can replace it with the canonical quadrature moment without leaking
 * that transitional representation into the FE shell.
 *
 * @tparam Space H(curl) primal polynomial space
 * @param space primal space on which the functional acts
 * @param tangent reference edge tangent
 * @param points functional support points
 * @param polynomial edge polynomial coefficients
 * @param functionals ordered destination dual block
 */
template<typename Space>
void
appendHCurlWeightedEdgeFunctional(
    Space const& space,
    typename node<typename Space::value_type>::type const& tangent,
    typename Space::points_type const& points,
    ublas::vector<typename Space::value_type> const& polynomial,
    std::vector<Functional<Space>>& functionals )
{
    functional::SecondDirectionalComponentPointEvaluation<Space> block(
        space, tangent, points, polynomial );
    std::copy( block.begin(), block.end(), std::back_inserter( functionals ) );
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HCURLFUNCTIONALS_HPP
