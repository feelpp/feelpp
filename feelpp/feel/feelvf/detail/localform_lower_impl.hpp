/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file localform_lower_impl.hpp
   \brief Lowering entry points for supported scalar local forms.
 */
#ifndef FEELPP_VF_DETAIL_LOCALFORM_LOWER_IMPL_HPP
#define FEELPP_VF_DETAIL_LOCALFORM_LOWER_IMPL_HPP 1

#include <feel/feelvf/detail/localform_expr.hpp>
#include <feel/feelvf/detail/localform_match_impl.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

template<typename ExprT>
auto lower_scalar_bilinear_localform( ExprT const& expr )
{
    using factors = scalar_localform_factors<ExprT>;
    static_assert( can_lower_scalar_bilinear_localform_v<ExprT>,
                   "lower_scalar_bilinear_localform only supports scalar coefficient * idt(u) * id(v) local forms" );
    using lowered_expr_type =
        LoweredScalarMassExpr<typename factors::test_space_type,
                              typename factors::trial_space_type,
                              typename factors::coefficient_expr_type>;
    return vf::expr( lowered_expr_type( factors::coefficient_expr( expr ), expr.polynomialOrder() ) );
}

template<typename ExprT>
auto lower_scalar_linear_localform( ExprT const& expr )
{
    using factors = scalar_localform_factors<ExprT>;
    static_assert( can_lower_scalar_linear_localform_v<ExprT>,
                   "lower_scalar_linear_localform only supports scalar coefficient * id(v) local forms" );
    using lowered_expr_type =
        LoweredScalarSourceExpr<typename factors::test_space_type,
                                typename factors::coefficient_expr_type>;
    return vf::expr( lowered_expr_type( factors::coefficient_expr( expr ), expr.polynomialOrder() ) );
}

} // namespace detail
} // namespace vf
} // namespace Feel

#endif
