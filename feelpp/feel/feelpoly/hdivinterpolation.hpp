/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP
#define FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP 1

#include <feel/feelpoly/basisinterpolation.hpp>
#include <feel/feelpoly/piolatransform.hpp>

#include <type_traits>

namespace Feel::detail
{

/**
 * @brief Apply facet-normal point DoFs to an expression.
 *
 * The normal provider is evaluated at every support point so that curved
 * geometric maps use the pointwise contravariant-Piola cofactor.  This kernel
 * is shared by RT and BDM, for both static and runtime orders.
 *
 * @tparam ExprType expression evaluated at functional support points
 * @tparam LocalVectorType mutable local coefficient vector
 * @tparam NormalType reusable normal-vector storage
 * @tparam NormalProviderType callable providing the physical Piola cofactor normal
 * @param expr expression to interpolate
 * @param localValues destination local coefficients
 * @param facetCount number of facets in the traversed entity
 * @param dofPerFacet number of point DoFs per facet
 * @param normal reusable normal workspace
 * @param normalProvider geometry policy called for every facet point
 */
template<typename ExprType, typename LocalVectorType, typename NormalType, typename NormalProviderType>
void
interpolateHDivFacetPointDofs( ExprType const& expr,
                               LocalVectorType& localValues,
                               uint16_type facetCount,
                               uint16_type dofPerFacet,
                               NormalType& normal,
                               NormalProviderType&& normalProvider )
{
    using shape_type = typename std::remove_cvref_t<ExprType>::shape;
    for ( uint16_type facet = 0; facet < facetCount; ++facet )
    {
        for ( uint16_type ordinal = 0; ordinal < dofPerFacet; ++ordinal )
        {
            const uint16_type point = static_cast<uint16_type>( facet*dofPerFacet + ordinal );
            normalProvider( facet, point, normal );
            for ( uint16_type component = 0; component < shape_type::M; ++component )
                localValues( point ) += expr.evalq( component, 0, point )*normal( component );
        }
    }
}

/** Apply the standard contravariant-Piola facet normals to point DoFs. */
template<typename ExprType, typename LocalVectorType, typename NormalType>
void
interpolateHDivFacetDofs( ExprType const& expr,
                          LocalVectorType& localValues,
                          uint16_type facetCount,
                          uint16_type dofPerFacet,
                          NormalType& normal )
{
    auto geometry = expr.geom();
    interpolateHDivFacetPointDofs(
        expr, localValues, facetCount, dofPerFacet, normal,
        [&geometry]( uint16_type facet, uint16_type point, auto& facetNormal )
        {
            const uint16_type elementFacet =
                geometry->faceId() == invalid_uint16_type_value
                    ? facet : geometry->faceId();
            geometry->piolaContravariantFacetNormal(
                elementFacet, point, facetNormal );
        } );
}

/**
 * Evaluate the reference-cell integrand used by an H(div) interior moment.
 *
 * For the contravariant Piola map
 *
 *     v = (1/J) K v_hat,
 *
 * the inverse pullback is v_hat = J K^{-1} v.  Consequently a reference
 * moment against q_hat is evaluated on the physical vector as
 *
 *     v_hat . q_hat = J v . (K^{-T} q_hat).
 *
 * Geometric expression contexts expose K^{-T} as B().  Lightweight
 * reference-element test expressions intentionally have no geometric map;
 * for those, the identity-map branch preserves the reference integrand.
 *
 * @tparam ExprType physical or reference expression type
 * @tparam PolynomialValuesType tabulated vector test-polynomial values
 * @param expr expression being interpolated
 * @param polynomialValues test-polynomial values at quadrature points
 * @param polynomialIndex selected vector polynomial
 * @param quadraturePoint quadrature-point index in polynomialValues
 * @param expressionPoint corresponding expression evaluation point
 * @param nComponents vector dimension
 * @return contravariant-Piola-correct interior-moment integrand
 */
template<typename ExprType, typename PolynomialValuesType>
auto
hdivInteriorMomentIntegrand( ExprType const& expr,
                             PolynomialValuesType const& polynomialValues,
                             int polynomialIndex,
                             int quadraturePoint,
                             int expressionPoint,
                             int nComponents )
{
    using value_type = std::decay_t<decltype( expr.evalq( 0, 0, expressionPoint ) )>;
    value_type result = 0;

    if constexpr ( requires {
                       expr.geom()->inverseJacobianTranspose( expressionPoint );
                       expr.geom()->integrationJacobian( expressionPoint );
                   } )
    {
        result = geometry::contravariantPiolaDualPairing(
            expr.geom()->inverseJacobianTranspose( expressionPoint ),
            expr.geom()->integrationJacobian( expressionPoint ),
            [&]( int component )
            { return expr.evalq( component, 0, expressionPoint ); },
            [&]( int component )
            {
                return polynomialValues(
                    nComponents * polynomialIndex + component, quadraturePoint );
            },
            nComponents, nComponents );
    }
    else
    {
        for ( int component = 0; component < nComponents; ++component )
            result += expr.evalq( component, 0, expressionPoint ) *
                      polynomialValues( nComponents * polynomialIndex + component, quadraturePoint );
    }

    return result;
}

/**
 * @brief Apply all H(div) interior moments for a prepared polynomial space.
 *
 * The caller supplies only the order-specific polynomial space and quadrature;
 * RT/BDM and static/runtime variants share evaluation, validation, Piola
 * pairing, and accumulation.
 */
template<typename ExprType, typename LocalVectorType,
         typename PolynomialSetType, typename QuadratureType>
void
interpolateHDivInteriorMomentDofs( ExprType const& expr,
                                   LocalVectorType& localValues,
                                   PolynomialSetType const& polynomialSet,
                                   QuadratureType const& quadrature,
                                   int firstInternalDof,
                                   uint16_type internalDofCount,
                                   uint16_type componentCount )
{
    CHECK_EQ( polynomialSet.polynomialDimension(), internalDofCount )
        << "Invalid H(div) interior interpolation dimension";
    if constexpr ( requires { expr.nPoints(); } )
    {
        CHECK( expr.nPoints() >= firstInternalDof + quadrature.nPoints() )
            << "H(div) interior-moment interpolation requires expression values at appended quadrature points: got "
            << expr.nPoints() << " points, need "
            << firstInternalDof + quadrature.nPoints();
    }

    auto const polynomialValues = polynomialSet.evaluate( quadrature.points() );
    for ( int polynomial = 0; polynomial < internalDofCount; ++polynomial )
    {
        const int dof = firstInternalDof + polynomial;
        for ( int point = 0; point < quadrature.nPoints(); ++point )
        {
            const int expressionPoint = firstInternalDof + point;
            const auto integrand = hdivInteriorMomentIntegrand(
                expr, polynomialValues, polynomial, point, expressionPoint,
                componentCount );
            localValues( dof ) += quadrature.weight( point ) * integrand;
        }
    }
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP
