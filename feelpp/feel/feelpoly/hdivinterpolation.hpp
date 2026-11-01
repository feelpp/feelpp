/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP
#define FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP 1

#include <feel/feelpoly/piolatransform.hpp>

#include <type_traits>

namespace Feel::detail
{

/**
 * @brief Present one FE test-basis row as an ordinary H(div) expression.
 *
 * This adapter lets interpolation-matrix assembly reuse the same normal and
 * interior-moment kernels as ordinary RT/BDM interpolation, with runtime
 * cardinalities taken from the destination Eigen matrix.
 */
template<typename BasisExpression>
class HDivBasisExpressionRow
{
public:
    using expression_type = std::remove_reference_t<BasisExpression>;
    using shape = typename expression_type::shape;

    /** Construct a view of one flattened test-basis row. */
    HDivBasisExpressionRow( expression_type& expression, uint16_type row )
        : M_expression( expression ), M_row( row )
    {}

    /** Return the mapped geometry context owned by the expression. */
    [[nodiscard]] decltype(auto) geom() const { return M_expression.geom(); }
    /** Return the number of functional support/quadrature samples. */
    [[nodiscard]] int nPoints() const { return M_expression.nPoints(); }
    /** Evaluate the selected test-basis row at one sample and component. */
    [[nodiscard]] decltype(auto) evalq( uint16_type c1,
                                        uint16_type c2,
                                        uint16_type q ) const
    {
        return M_expression.evaliq( M_row, c1, c2, q );
    }

private:
    expression_type& M_expression;
    uint16_type M_row;
};

/** Construct a zero-copy H(div) test-basis row expression. */
template<typename BasisExpression>
[[nodiscard]] auto
hdivBasisExpressionRow( BasisExpression& expression, uint16_type row )
{
    return HDivBasisExpressionRow<BasisExpression>( expression, row );
}

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

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP
