/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
   \file sb9_stabilization.hpp
   \brief SB9 stabilization shell operators
 */
#ifndef FEELPP_VF_SB9_STABILIZATION_HPP
#define FEELPP_VF_SB9_STABILIZATION_HPP 1

#include <feel/feelvf/sb9_common.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{
enum class SB9StabilizationKind
{
    Bs1,
    Bs2,
    Bs3,
    Bs4
};

template <typename GeometryDataType>
class SB9StabilizationKernelCache : public SB9KernelBase<GeometryDataType>
{
public:
    using base_type = SB9KernelBase<GeometryDataType>;
    static constexpr uint16_type node_count = base_type::node_count;
    static constexpr uint16_type component_count = base_type::component_count;

    using base_type::base_type;

    template <SB9StabilizationKind Kind, typename VectorType>
    void fillStabilizationCoefficients( VectorType& coeff, uint16_type node, uint16_type component ) const
    {
        static_assert( Kind == SB9StabilizationKind::Bs1 || Kind == SB9StabilizationKind::Bs2 ||
                       Kind == SB9StabilizationKind::Bs3 || Kind == SB9StabilizationKind::Bs4,
                       "unsupported SB9 stabilization kind" );

        auto const& frame = this->M_frame[component];

        if constexpr ( Kind == SB9StabilizationKind::Bs1 )
        {
            coeff( 0 ) = frame[2] * this->M_data.vgamma( node, 0 );
        }
        else if constexpr ( Kind == SB9StabilizationKind::Bs2 )
        {
            coeff( 0 ) = frame[2] * this->M_data.vgamma( node, 1 );
        }
        else if constexpr ( Kind == SB9StabilizationKind::Bs3 )
        {
            coeff( 0 ) = frame[0] * this->M_data.vgamma( node, 2 );
            coeff( 1 ) = frame[1] * this->M_data.vgamma( node, 2 );
        }
        else
        {
            coeff( 0 ) = frame[0] * this->M_data.vgamma( node, 3 );
            coeff( 1 ) = frame[1] * this->M_data.vgamma( node, 3 );
            coeff( 2 ) = frame[2] * this->M_data.vgamma( node, 3 );
        }
    }
};
} // namespace detail

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs1( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                       proxy_type::role,
                                                       detail::SB9StabilizationKernelCache,
                                                       detail::SB9StabilizationKind::Bs1,
                                                       1>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs2( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                       proxy_type::role,
                                                       detail::SB9StabilizationKernelCache,
                                                       detail::SB9StabilizationKind::Bs2,
                                                       1>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs3( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                       proxy_type::role,
                                                       detail::SB9StabilizationKernelCache,
                                                       detail::SB9StabilizationKind::Bs3,
                                                       2>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9Bs4( ProxyType const& proxy )
{
    using proxy_type = detail::basis_proxy_type_t<ProxyType>;
    using expr_type = detail::SB9StabilizationOperator<typename proxy_type::element_type,
                                                       proxy_type::role,
                                                       detail::SB9StabilizationKernelCache,
                                                       detail::SB9StabilizationKind::Bs4,
                                                       3>;
    return Expr<expr_type>( expr_type( proxy.element() ) );
}
} // namespace vf
} // namespace Feel

#endif
