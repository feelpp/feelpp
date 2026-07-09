/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
    SPDX-FileContributor: Hanna Chetouane

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
   \file sb9_pinching.hpp
   \brief SB9 pinching shell kinematic operators
 */
#ifndef FEELPP_VF_SB9_PINCHING_HPP
#define FEELPP_VF_SB9_PINCHING_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
/**
 * \brief Selects the SB9 pinching coefficient family.
 *
 * The values are used as compile-time tags by \ref SB9PinchingKernelCache and
 * \ref SB9VectorOperator.
 */
enum class SB9PinchingKind
{
    /// Constant-throught-thickness pinching coefficient vector, usually denoted Bpc.
    Bpc,
    /// Linear-through-thickness pinching coefficient vector, usually denoted Bpz.
    Bpz
};

/**
 * \brief Element-local cache for SB9 pinching coefficients.
 *
 * The cache reuses \ref SB9KernelBase for frame and Jacobian data and computes
 * the pinching coefficients. The resulting coefficients are later projected into
 * Mandel symmetric storage by \ref fillVectorCoefficients.
 *
 * \tparam GeometryDataType Shell cell geometry cache type provided by
 *         ShellCellGeometryTensorBase.
 */
template <typename GeometryDataType>
class SB9PinchingKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    /// Shared SB9 geometry cache base.
    using base_type = SB9KernelBase<GeometryDataType>;
    /// Scalar type used by the geometry and generated coefficients.
    using value_type = typename base_type::value_type;
    /// Number of geometric nodes in the current SB9 Q1 hexahedral element.
    static constexpr uint16_type node_count = base_type::node_count;
    /// Number of displacement components handled by SB9 pinching operators.
    static constexpr uint16_type component_count = base_type::component_count;

    /**
     * \brief Build pinching coefficient data for one element.
     *
     * This constructor reuses the element geometry data to precompute the
     * linear-through-thickness pinching coefficients for each node required by the SB9 
     * pinching operator.
     *
     * \param data Element-local shell geometry data.
     */
    explicit SB9PinchingKernelCache( GeometryDataType const& data )
        :
        base_type( data )
    {
        for ( uint16_type node = 0; node < node_count; ++node )
        {
            M_pinching[node] = this->M_data.vgamma( node, 0 )*this->M_invJ0( 2, 1 ) +
                               this->M_data.vgamma( node, 1 )*this->M_invJ0( 2, 0 );
        }
    }

    /**
     * \brief Fill one SB9 pinching coefficient vector.
     *
     * \tparam Kind Compile-time selector, either \ref SB9PinchingKind::Bpc or
     *         \ref SB9PinchingKind::Bpz.
     * \tparam VectorType Eigen-compatible coefficient vector type.
     * \param coeff Output coefficient vector in Feel++ symmetric storage order.
     * \param node Local geometric node index.
     * \param component Displacement component index.
     */
    template <SB9PinchingKind Kind, typename VectorType>
    void fillVectorCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9PinchingKind::Bpc || Kind == SB9PinchingKind::Bpz,
                       "unsupported SB9 pinching vector kind" );

        if constexpr ( Kind == SB9PinchingKind::Bpc )
            this->fillPinchingCoefficients( coeff, component, this->M_data.bz( node ) );
        else
            this->fillPinchingCoefficients( coeff, component, M_pinching[node] );            
    }

private:
    /// Precomputed the linear-through-thickness pinching coefficients per node.
    std::array<value_type, node_count> M_pinching{};
};
} // namespace detail

/**
 * \brief Build the constant-through-thickness SB9 pinching coefficient expression.
 *
 * The returned expression evaluates the `Bpc` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bpc` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bpc( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9PinchingKernelCache,
                                                detail::SB9PinchingKind::Bpc>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the linear-through-thickness SB9 pinching coefficient expression.
 *
 * The returned expression evaluates the `Bpz` contribution for a trial or test
 * basis proxy created with `trial(Vh, ...)` or `test(Vh, ...)`.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \param proxy Trial or test basis proxy.
 * \return Feel++ expression containing the SB9 `Bpz` coefficients.
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bpz( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9VectorOperator<typename proxy_type::element_type,
                                                proxy_type::role,
                                                detail::SB9PinchingKernelCache,
                                                detail::SB9PinchingKind::Bpz>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

/**
 * \brief Build the SB9 pinching Mandel strain expression associated with the Q1
 * displacement field.
 *
 * This helper combines the pinching coefficient expressions as
 * `Bpc + bpzScale * zeta * Bpz` and returns a six-component Mandel vector,
 * with all other components set to zero.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ZetaExprT Feel++ expression type used for the through-thickness
 *         coordinate, typically `zeta()`.
 * \tparam ScaleExprT Feel++ expression type used to scale the `zeta * Bpz`
 *         contribution.
 * \param proxy Trial or test basis proxy.
 * \param zetaExpr Through-thickness coordinate or scaling expression.
 * \param bpzScaleExpr Scale applied to the linear-through-thickness `Bpz`
 *        contribution.
 * \return Feel++ Mandel-vector expression for the SB9 pinching strain.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT, typename ScaleExprT>
[[nodiscard]] inline auto
sb9Pinching( ProxyType const& proxy, ZetaExprT const& zetaExpr, ScaleExprT const& bpzScaleExpr )
{
    auto bpc = sb9Bpc( proxy );
    auto bpz = sb9Bpz( proxy );

    return mandel_component<3,2,2>( component<2,0>( bpc ) + bpzScaleExpr * zetaExpr * component<2,0>( bpz ) );
}

/**
 * \brief Build the SB9 pinching Mandel strain expression associated with the Q1
 * displacement field.
 *
 * This helper combines the pinching coefficient expressions as `Bpc + zeta * Bpz`
 * and returns a six-component Mandel vector, with all other components set to zero.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ZetaExprT Feel++ expression type used for the through-thickness
 *         coordinate, typically `zeta()`.
 * \param proxy Trial or test basis proxy.
 * \param zetaExpr Through-thickness coordinate or scaling expression.
 * \return Feel++ Mandel-vector expression for the SB9 pinching strain.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT>
[[nodiscard]] inline auto
sb9Pinching( ProxyType const& proxy, ZetaExprT const& zetaExpr )
{
    return sb9Pinching( proxy, zetaExpr, cst( 1.0 ) );
}

/**
 * \brief Build the SB9 pinching Mandel strain expression using the default
 * through-thickness coordinate expression \ref zeta().
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Pinching( ProxyType const& proxy )
{
    return sb9Pinching( proxy, zeta() );
}

/**
 * \brief Build the SB9 pinching Mandel strain expression associated with the
 * scalar space carrying the 25th SB9 degree of freedom.
 *
 * This helper builds the scalar pinching contribution as `-4 * zeta / h0` in
 * the pinching Mandel entry and returns a six-component Mandel vector with all
 * other components set to zero.
 *
 * \tparam ProxyType Feel++ trial/test basis proxy type.
 * \tparam ZetaExprT Feel++ expression type used for the through-thickness
 *         coordinate, typically `zeta()`.
 * \param proxy Trial or test basis proxy.
 * \param zetaExpr Through-thickness coordinate or scaling expression.
 * \return Feel++ Mandel-vector expression for the SB9 pinching strain.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT>
[[nodiscard]] inline auto
sb9PinchingW9( ProxyType const& proxy, ZetaExprT const& zetaExpr )
{
    return mandel_component<3,2,2>( ( -4.0 * zetaExpr / shellThickness() ) * id( proxy ) );
}

/**
 * \brief Build the SB9 pinching Mandel strain expression using the default
 * through-thickness coordinate expression \ref zeta().
 */
template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9PinchingW9( ProxyType const& proxy )
{
    return sb9PinchingW9( proxy, zeta() );
}
} // namespace vf
} // namespace Feel

#endif
