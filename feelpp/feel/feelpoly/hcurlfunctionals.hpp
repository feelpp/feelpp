/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/** @file hcurlfunctionals.hpp @brief H(curl) family dual-functional builders. */

#ifndef FEELPP_FEELPOLY_HCURLFUNCTIONALS_HPP
#define FEELPP_FEELPOLY_HCURLFUNCTIONALS_HPP 1

#include <feel/feelpoly/tracefunctionals.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>
#include <feel/feelpoly/hcurlmomentpolicy.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
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
template<typename Space, typename ReferenceConvexType, typename EdgePointsProvider>
void
appendHCurlEdgeTangentialPointFunctionals(
    Space const& space,
    ReferenceConvexType const& referenceConvex,
    uint16_type dofPerEdge,
    typename Space::points_type& elementPoints,
    std::vector<typename Space::points_type>& facetPoints,
    std::vector<Functional<Space>>& functionals,
    EdgePointsProvider&& edgePointsProvider )
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
        points_type edgePoints( std::invoke( edgePointsProvider, edge ) );

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

/** @brief Append edge functionals using the reference convex point provider. */
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
    appendHCurlEdgeTangentialPointFunctionals(
        space, referenceConvex, dofPerEdge, elementPoints, facetPoints,
        functionals,
        [&referenceConvex]( uint16_type edge )
        { return referenceConvex.makePoints( 1, edge ); } );
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

/**
 * Append canonical tangential polynomial moments on every triangular face of
 * a reference tetrahedron.
 *
 * For each face, two edge vectors from its first vertex define a covariant
 * tangential frame and affine coordinates `(a,b)`. Moments are ordered by
 * face, tangent, total monomial degree, then the exponent of `b`.
 *
 * @param space tetrahedral H(curl) primal space
 * @param referenceConvex reference tetrahedron
 * @param polynomialDegree degree of the scalar face test space
 * @param functionals ordered destination dual block
 */
template<typename Space, typename ReferenceConvexType>
void
appendHCurlTriangularFaceMoments( Space const& space,
                                  ReferenceConvexType const& referenceConvex,
                                  int polynomialDegree,
                                  std::vector<Functional<Space>>& functionals )
{
    static_assert( Space::nDim == 3 );
    if ( polynomialDegree < 0 )
        return;

    using value_type = typename Space::value_type;
    // Runtime quadrature is required when Space::nOrder is Dynamic.  The
    // requested degree is the face test degree k-1, hence 2*(k+2) preserves
    // the quadrature envelope used by the static NED1 implementation.
    IMGeneral<3, value_type, Simplex> quadrature(
        hcurlTriangularFaceMomentQuadratureOrder( polynomialDegree ) );
    for ( uint16_type face = 0; face < 4; ++face )
    {
        auto const vertices = referenceConvex.faceVertices( face );
        Eigen::Matrix<value_type,3,1> origin, tangent0, tangent1;
        for ( int c = 0; c < 3; ++c )
        {
            origin( c ) = vertices( c, 0 );
            tangent0( c ) = vertices( c, 1 ) - vertices( c, 0 );
            tangent1( c ) = vertices( c, 2 ) - vertices( c, 0 );
        }
        const value_type g00 = tangent0.dot( tangent0 );
        const value_type g01 = tangent0.dot( tangent1 );
        const value_type g11 = tangent1.dot( tangent1 );
        const value_type metricDeterminant = g00*g11-g01*g01;
        const value_type inverseDeterminant = value_type( 1 )/metricDeterminant;
        const value_type surfaceJacobian = std::sqrt( metricDeterminant );
        auto const& points = quadrature.fpoints( face );
        auto const basisAtPoints = functional::detail::basisEvaluateAtPoints( space, points );

        for ( auto const& tangent : { tangent0, tangent1 } )
        {
            for ( int degree = 0; degree <= polynomialDegree; ++degree )
            {
                for ( int exponentB = 0; exponentB <= degree; ++exponentB )
                {
                    const int exponentA = degree-exponentB;
                    typename Space::matrix_type coefficients( Space::nComponents,
                                                               basisAtPoints.size1() );
                    coefficients.clear();
                    for ( uint16_type q = 0; q < quadrature.nPointsOnFace( face ); ++q )
                    {
                        Eigen::Matrix<value_type,3,1> delta;
                        for ( int c = 0; c < 3; ++c )
                            delta( c ) = points( c, q )-origin( c );
                        const value_type rhs0 = tangent0.dot( delta );
                        const value_type rhs1 = tangent1.dot( delta );
                        const value_type a = ( g11*rhs0-g01*rhs1 )*inverseDeterminant;
                        const value_type b = ( g00*rhs1-g01*rhs0 )*inverseDeterminant;
                        // Use the intrinsic reference-triangle measure.  The
                        // embedded tetrahedron faces have different surface
                        // Jacobians, which must not enter a shared trace DoF.
                        const value_type weight = quadrature.weight( face, q )/surfaceJacobian*
                                                  std::pow( a, exponentA )*std::pow( b, exponentB );
                        for ( uint16_type component = 0; component < 3; ++component )
                            for ( uint16_type basis = 0; basis < basisAtPoints.size1(); ++basis )
                                coefficients( component, basis ) +=
                                    weight*tangent( component )*basisAtPoints( basis, q );
                    }
                    functionals.emplace_back( space, coefficients );
                }
            }
        }
    }
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HCURLFUNCTIONALS_HPP
