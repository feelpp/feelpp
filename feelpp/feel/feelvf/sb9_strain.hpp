/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

/**
   \file sb9_strain.hpp
   \brief SB9 shell strain composition helpers
 */
#ifndef FEELPP_VF_SB9_STRAIN_HPP
#define FEELPP_VF_SB9_STRAIN_HPP 1

#include <feel/feelvf/sb9_bending.hpp>
#include <feel/feelvf/sb9_pinching.hpp>
#include <feel/feelvf/sb9_shear.hpp>

namespace Feel
{
namespace vf
{
/**
 * \brief Compose the complete displacement-driven SB9 strain in Mandel storage.
 *
 * \param proxy Trial or test displacement basis proxy.
 * \param zetaExpr Through-thickness coordinate expression.
 * \param shearWeight Reissner transverse-shear weight.
 * \param pinchingBpzScale Scale applied to the linear pinching contribution.
 * \return SB9 strain in Feel++ storage order `(00,01,02,11,12,22)`.
 */
template <detail::BasisProxyType ProxyType,
          typename ZetaExprT,
          typename ShearWeightExprT,
          typename PinchingScaleExprT>
[[nodiscard]] inline auto
sb9ShellStrain( ProxyType const& proxy,
                ZetaExprT const& zetaExpr,
                ShearWeightExprT const& shearWeight,
                PinchingScaleExprT const& pinchingBpzScale )
{
    auto membraneBending = sb9MembraneBending( proxy, zetaExpr );
    auto pinching = sb9Pinching( proxy, zetaExpr, pinchingBpzScale );
    auto shear = sb9Shear( proxy, shearWeight );

    // Rebuild the final symmetric-storage vector explicitly in Feel's
    // internal storage order (00,01,02,11,12,22).
    return vec( component<0, 0>( membraneBending ),
                component<1, 0>( membraneBending ),
                component<2, 0>( shear ),
                component<3, 0>( membraneBending ),
                component<4, 0>( shear ),
                component<5, 0>( pinching ) );
}

/**
 * \brief Compose the complete SB9 strain with the full `zeta*Bpz` contribution.
 */
template <detail::BasisProxyType ProxyType, typename ZetaExprT, typename ShearWeightExprT>
[[nodiscard]] inline auto
sb9ShellStrain( ProxyType const& proxy,
                ZetaExprT const& zetaExpr,
                ShearWeightExprT const& shearWeight )
{
    return sb9ShellStrain( proxy, zetaExpr, shearWeight, cst( 1.0 ) );
}

template <detail::BasisProxyType ProxyType>
[[nodiscard]] inline auto
sb9W9( ProxyType const& proxy )
{
    return mandel_component<3, 2, 2>( ( -4.0*zeta()/shellThickness() )*id( proxy ) );
}

template <detail::BasisProxyType ProxyType, typename ScaleExprT>
[[nodiscard]] inline auto
sb9W9( ProxyType const& proxy, ScaleExprT const& scaleExpr )
{
    return scale_symm_storage( scaleExpr, sb9W9( proxy ) );
}
} // namespace vf
} // namespace Feel

#endif
