/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file sb9_bending.hpp
   \brief SB9 membrane and bending shell kinematic operators
 */
#ifndef FEELPP_VF_SB9_BENDING_HPP
#define FEELPP_VF_SB9_BENDING_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Selects the SB9 membrane or bending coefficient family.
 *
 * The values are used as compile-time tags by \ref SB9BendingKernelCache and
 * \ref SB9VectorOperator.
 */
enum class SB9BendingKind
{
    /// Mid-surface membrane coefficient matrix, usually denoted Bm0.
    Bm0,
    /// Linear-through-thickness bending coefficient matrix, usually denoted Bb0.
    Bb0
};

/**
 * \brief Element-local cache for SB9 membrane and bending coefficients.
 *
 * The cache reuses \ref SB9KernelBase for frame and Jacobian data and adds the
 * bending derivative coefficients derived from the shell geometry `vgamma`
 * modes. The coefficients are later projected into Mandel symmetric storage by
 * \ref fillVectorCoefficients.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9BendingKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    /// Shared SB9 geometry cache base.
    using base_type = SB9KernelBase<GeometryDataType>;
    /// Scalar type used by the geometry and generated coefficients.
    using value_type = typename base_type::value_type;
    /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = base_type::node_count;
    /// Number of displacement components handled by SB9 bending operators.
    static constexpr uint16_type component_count = base_type::component_count;

    /**
     * \brief Build membrane/bending coefficient data for one element.
     *
     * The membrane part is available directly from the geometry cache through
     * `bx` and `by`; this constructor precomputes the two bending derivative
     * coefficients for each node.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9BendingKernelCache( GeometryDataType const& data )
        :
        base_type( data )
    {
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            M_bending[node][0] = this->M_data.vgamma( node, 1 )*this->M_invJ0( 0, 0 );
            M_bending[node][1] = this->M_data.vgamma( node, 0 )*this->M_invJ0( 1, 1 ) +
                                 this->M_data.vgamma( node, 1 )*this->M_invJ0( 1, 0 );
        }
    }

    /**
     * \brief Fill one SB9 bending-family coefficient vector.
     *
     * \tparam Kind Compile-time selector, either \ref SB9BendingKind::Bm0 or
     *         \ref SB9BendingKind::Bb0.
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in Feel++ symmetric storage order.
     * \param node Local geometric node index.
     * \param component Displacement component index.
     */
    template <SB9BendingKind Kind, typename VectorType>
    void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9BendingKind::Bm0 || Kind == SB9BendingKind::Bb0,
                       "unsupported SB9 bending vector kind" );

        if constexpr ( Kind == SB9BendingKind::Bm0 )
            this->fillMembraneCoefficients( coeff, component, this->M_data.bx( node ), this->M_data.by( node ) );
            // this->fillMembraneCoefficients( coeff, component, value_type( 12.0 ), value_type( 12.0 ) );
        else
            this->fillMembraneCoefficients( coeff, component, M_bending[node][0], M_bending[node][1] );
            // this->fillMembraneCoefficients( coeff, component, value_type( 12.0 ), value_type( 12.0 ) );
    }

private:
    /// Precomputed bending derivative coefficients per node and in-plane axis.
    std::array<std::array<value_type, 2>, node_count> M_bending{};
};
} // namespace detail

/**
 * \brief Build the SB9 mid-surface membrane coefficient expression.
 *
 * The returned expression evaluates the `Bm0` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bm0` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bm0( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9BendingKernelCache,
                                                detail::SB9BendingKind::Bm0>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 bending coefficient expression.
 *
 * The returned expression evaluates the `Bb0` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bb0` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bb0( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9BendingKernelCache,
                                                detail::SB9BendingKind::Bb0>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 membrane-plus-bending Mandel strain expression.
 *
 * This helper combines the membrane and bending coefficient expressions as
 * `Bm0 + zeta * Bb0` in the in-plane symmetric-storage entries and returns a
 * six-component Mandel vector with out-of-plane components set to zero.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ZetaExprT Feel++ expression type used for the through-thickness
 *         coordinate, typically `zeta()`.
 * \param proxy Trial or test basis proxy.
 * \param zetaExpr Through-thickness coordinate or scaling expression.
 * \return Feel++ Mandel-vector expression for the SB9 membrane/bending strain.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT>
[[nodiscard]] inline auto
sb9MembraneBending( ProxyType const& proxy, ZetaExprT const& zetaExpr )
{
    auto bm0 = sb9Bm0( proxy );
    auto bb0 = sb9Bb0( proxy );

    return mandel_vec<3>( component<0, 0>( bm0 ) + zetaExpr * component<0, 0>( bb0 ),
                          component<1, 0>( bm0 ) + zetaExpr * component<1, 0>( bb0 ),
                          cst( 0.0 ),
                          component<3, 0>( bm0 ) + zetaExpr * component<3, 0>( bb0 ),
                          cst( 0.0 ),
                          cst( 0.0 ) );
}
} // namespace vf
} // namespace Feel

#endif
