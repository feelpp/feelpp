/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#ifndef FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP
#define FEELPP_FEELPOLY_HDIVINTERPOLATION_HPP 1

#include <type_traits>

namespace Feel::detail
{

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

    if constexpr ( requires { expr.geom()->B( expressionPoint ); expr.geom()->J( expressionPoint ); } )
    {
        auto const& inverseJacobianTranspose = expr.geom()->B( expressionPoint );
        for ( int physicalComponent = 0; physicalComponent < nComponents; ++physicalComponent )
        {
            value_type transformedTest = 0;
            for ( int referenceComponent = 0; referenceComponent < nComponents; ++referenceComponent )
                transformedTest += inverseJacobianTranspose( physicalComponent, referenceComponent ) *
                                   polynomialValues( nComponents * polynomialIndex + referenceComponent,
                                                     quadraturePoint );
            result += expr.evalq( physicalComponent, 0, expressionPoint ) * transformedTest;
        }
        result *= expr.geom()->J( expressionPoint );
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
