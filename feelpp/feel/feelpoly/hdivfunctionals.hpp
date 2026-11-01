/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELPOLY_HDIVFUNCTIONALS_HPP
#define FEELPP_FEELPOLY_HDIVFUNCTIONALS_HPP 1

#include <feel/feelpoly/tracefunctionals.hpp>
#include <feel/feelpoly/pointsetquadrature.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

namespace Feel::detail
{

/**
 * @brief Compute the measure of a reference-simplex facet.
 * @tparam ReferenceConvexType triangle or tetrahedron reference type
 * @param referenceConvex reference cell
 * @param facet local facet index
 * @return edge length in 2D or triangle area in 3D
 */
template<typename ReferenceConvexType>
auto
referenceSimplexFacetMeasure( ReferenceConvexType const& referenceConvex, uint16_type facet )
{
    using value_type = typename ReferenceConvexType::value_type;
    static constexpr uint16_type nDim = ReferenceConvexType::nDim;
    static_assert( nDim == 2 || nDim == 3,
                   "H(div) simplex facet functionals currently support triangles and tetrahedra." );

    auto const vertices = referenceConvex.faceVertices( facet );
    if constexpr ( nDim == 2 )
    {
        value_type squaredLength = 0;
        for ( uint16_type c = 0; c < nDim; ++c )
        {
            const value_type delta = vertices( c, 1 ) - vertices( c, 0 );
            squaredLength += delta * delta;
        }
        return std::sqrt( squaredLength );
    }
    else
    {
        const value_type ax = vertices( 0, 1 ) - vertices( 0, 0 );
        const value_type ay = vertices( 1, 1 ) - vertices( 1, 0 );
        const value_type az = vertices( 2, 1 ) - vertices( 2, 0 );
        const value_type bx = vertices( 0, 2 ) - vertices( 0, 0 );
        const value_type by = vertices( 1, 2 ) - vertices( 1, 0 );
        const value_type bz = vertices( 2, 2 ) - vertices( 2, 0 );
        const value_type cx = ay*bz - az*by;
        const value_type cy = az*bx - ax*bz;
        const value_type cz = ax*by - ay*bx;
        return value_type( 0.5 ) * std::sqrt( cx*cx + cy*cy + cz*cz );
    }
}

/**
 * @brief Append nodal normal-flux functionals on all simplex facets.
 * @param space H(div) primal polynomial space
 * @param referenceConvex reference cell providing scaled facet normals
 * @param facetPoints support points grouped by facet
 * @param functionals ordered destination dual block
 */
template<typename Space, typename ReferenceConvexType, typename FacetPointsType>
void
appendHDivFacetNormalPointFunctionals( Space const& space,
                                       ReferenceConvexType const& referenceConvex,
                                       FacetPointsType const& facetPoints,
                                       std::vector<Functional<Space>>& functionals )
{
    using value_type = typename Space::value_type;
    using node_type = typename ReferenceConvexType::node_type;

    for ( int facet = referenceConvex.entityRange( Space::nDim-1 ).begin();
          facet < referenceConvex.entityRange( Space::nDim-1 ).end();
          ++facet )
    {
        node_type direction( Space::nDim );
        em_node_type<value_type> mappedDirection( direction.data().begin(), direction.size() );
        mappedDirection = referenceConvex.normal( facet ) *
                          referenceSimplexFacetMeasure( referenceConvex, facet );
        functional::appendDirectionalPointFunctionals( space, direction,
                                                       facetPoints[facet], functionals );
    }
}

/**
 * @brief Build a reference-cell L2 moment against a vector polynomial.
 * @param space H(div) primal polynomial space
 * @param polynomial vector test polynomial
 * @param quadratureOrder polynomial order used to select exact quadrature
 * @return mathematical functional represented in the current primal basis
 */
template<typename Space, typename Polynomial>
[[nodiscard]] Functional<Space>
makeHDivIntegralMomentFunctional( Space const& space,
                                  Polynomial const& polynomial,
                                  uint16_type quadratureOrder )
{
    using value_type = typename Space::value_type;
    static constexpr uint16_type nDim = Space::nDim;

    IMGeneral<nDim, value_type, Simplex> im( 2*quadratureOrder );
    auto const basisAtQuadPts = functional::detail::basisEvaluateAtPoints( space, im.points() );
    auto const polynomialAtQuadPts = polynomial.evaluate( im.points() );

    typename Space::matrix_type coeff( Space::nComponents, basisAtQuadPts.size1() );
    coeff.clear();
    for ( uint16_type component = 0; component < Space::nComponents; ++component )
    {
        for ( uint16_type basis = 0; basis < basisAtQuadPts.size1(); ++basis )
        {
            value_type value = 0;
            for ( uint16_type q = 0; q < im.nPoints(); ++q )
                value += im.weight( q ) * polynomialAtQuadPts( component, q ) * basisAtQuadPts( basis, q );
            coeff( component, basis ) = value;
        }
    }
    return Functional<Space>( space, coeff );
}

/**
 * Build an exact L2 moment when the test polynomial and primal space share the
 * same orthonormal scalar expansion basis (the BDM interior-moment case).
 * @param space H(div) primal polynomial space
 * @param polynomial vector test polynomial in the shared expansion basis
 * @return mathematical functional represented in the current primal basis
 */
template<typename Space, typename Polynomial>
[[nodiscard]] Functional<Space>
makeHDivOrthonormalIntegralMomentFunctional( Space const& space,
                                             Polynomial const& polynomial )
{
    typename Space::matrix_type coefficients( Space::nComponents, space.coeff().size2() );
    coefficients.clear();
    auto const& polynomialCoefficients = polynomial.coeff();
    CHECK( polynomialCoefficients.size1() == Space::nComponents )
        << "invalid H(div) moment polynomial component count";
    CHECK( polynomialCoefficients.size2() <= coefficients.size2() )
        << "invalid H(div) moment polynomial basis dimension";
    ublas::subrange( coefficients, 0, Space::nComponents,
                     0, polynomialCoefficients.size2() ) = polynomialCoefficients;
    return Functional<Space>( space, coefficients );
}

/**
 * @brief Build the complete ordered dual of an H(div) simplex element.
 *
 * RT and BDM differ only in their interior test space and moment functional.
 * Those policies are supplied by the caller; facet sampling, point layout,
 * cardinality checks, and functional ordering remain one mathematical path
 * for static and runtime orders.
 */
template<typename Space, typename ReferenceConvex,
         typename FacetPointsProvider, typename InteriorSpaceProvider,
         typename InteriorFunctionalProvider>
[[nodiscard]] std::vector<Functional<Space>>
buildHDivSimplexFunctionals(
    Space const& primal,
    ReferenceConvex const& referenceConvex,
    uint16_type facetDof,
    uint16_type interiorDof,
    uint16_type totalDof,
    uint16_type interiorQuadratureOrder,
    typename Space::points_type& elementPoints,
    std::vector<typename Space::points_type>& facetPoints,
    FacetPointsProvider&& facetPointsProvider,
    InteriorSpaceProvider&& interiorSpaceProvider,
    InteriorFunctionalProvider&& interiorFunctionalProvider )
{
    static constexpr uint16_type Dim = Space::nDim;
    using value_type = typename Space::value_type;
    using convex_type = typename ReferenceConvex::super;

    facetPoints.resize( convex_type::numTopologicalFaces );
    uint32_type pointOffset = 0;
    for ( int facet = referenceConvex.entityRange( Dim-1 ).begin();
          facet < referenceConvex.entityRange( Dim-1 ).end(); ++facet )
    {
        auto points = facetPointsProvider( static_cast<uint16_type>( facet ) );
        CHECK_EQ( points.size2(), facetDof )
            << "invalid H(div) facet point count on facet " << facet;
        facetPoints[facet] = points;
        const uint32_type nextPointOffset = pointOffset+points.size2();
        if ( elementPoints.size2() < nextPointOffset )
            elementPoints.resize( Dim, nextPointOffset, true );
        if ( points.size2() )
            ublas::subrange( elementPoints, 0, Dim,
                             pointOffset, nextPointOffset ) = points;
        pointOffset = nextPointOffset;
    }

    std::vector<Functional<Space>> functionals;
    functionals.reserve( totalDof );
    appendHDivFacetNormalPointFunctionals(
        primal, referenceConvex, facetPoints, functionals );

    if ( interiorDof > 0 )
    {
        auto interiorSpace = interiorSpaceProvider();
        CHECK_EQ( interiorSpace.polynomialDimension(), interiorDof )
            << "invalid H(div) interior moment count";
        IMGeneral<Dim, value_type, Simplex> quadrature(
            interiorQuadratureOrder );
        const uint32_type nextPointOffset = pointOffset+quadrature.nPoints();
        if ( elementPoints.size2() < nextPointOffset )
            elementPoints.resize( Dim, nextPointOffset, true );
        ublas::subrange( elementPoints, 0, Dim,
                         pointOffset, nextPointOffset ) = quadrature.points();
        for ( int i = 0; i < interiorSpace.polynomialDimension(); ++i )
            functionals.push_back( interiorFunctionalProvider(
                interiorSpace.polynomial( i ) ) );
    }

    CHECK_EQ( functionals.size(), totalDof )
        << "invalid H(div) functional count";
    return functionals;
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HDIVFUNCTIONALS_HPP
