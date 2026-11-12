/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELPOLY_BASISINTERPOLATION_HPP
#define FEELPP_FEELPOLY_BASISINTERPOLATION_HPP 1

#include <type_traits>

#include <feel/feelcore/feeltypes.hpp>

namespace Feel::detail
{

/**
 * @brief Present one FE test-basis row as an ordinary interpolation expression.
 *
 * Dual interpolation kernels operate on `evalq`, whereas interpolation-matrix
 * assembly evaluates every domain test basis through `evaliq`. This zero-copy
 * adapter is independent of the H(div) or H(curl) mapping policy.
 */
template<typename BasisExpression>
class BasisExpressionRow
{
public:
    using expression_type = std::remove_reference_t<BasisExpression>;
    using shape = typename expression_type::shape;

    BasisExpressionRow( expression_type& expression, uint16_type row )
        : M_expression( expression ), M_row( row )
    {}

    [[nodiscard]] decltype(auto) geom() const { return M_expression.geom(); }
    [[nodiscard]] int nPoints() const { return M_expression.nPoints(); }
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

/** Construct a zero-copy test-basis row expression. */
template<typename BasisExpression>
[[nodiscard]] auto
basisExpressionRow( BasisExpression& expression, uint16_type row )
{
    return BasisExpressionRow<BasisExpression>( expression, row );
}

/**
 * @brief Apply an FE's interpolation kernels to every domain basis row.
 *
 * Static and runtime FE wrappers use the same matrix materialization. The FE
 * remains responsible only for allocating one local vector and applying its
 * family-specific interpolation kernels.
 */
template<typename FiniteElementType, typename BasisExpression, typename LocalMatrix>
void
interpolateBasisFunctionRows( FiniteElementType const& finiteElement,
                              BasisExpression& expression,
                              LocalMatrix& localMatrix )
{
    localMatrix.setZero();
    CHECK_EQ( localMatrix.cols(), finiteElement.localDof() );
    for ( uint16_type row = 0;
          row < static_cast<uint16_type>( localMatrix.rows() ); ++row )
    {
        auto rowExpression = basisExpressionRow( expression, row );
        auto rowValues = finiteElement.localInterpolant();
        finiteElement.interpolate( rowExpression, rowValues );
        localMatrix.row( row ) = rowValues.transpose();
    }
}

} // namespace Feel::detail

#endif // FEELPP_FEELPOLY_BASISINTERPOLATION_HPP
