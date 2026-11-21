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

/**
 * @brief Lower a scalar isotropic H1 stiffness expression.
 * @param expr expression matching coefficient times gradient inner product
 * @return expression sharing the canonical local-tensor implementation
 */
template<typename ExprT>
auto lower_scalar_h1_stiffness_localform( ExprT const& expr )
{
    using factors = scalar_h1_stiffness_factors<ExprT>;
    static_assert(
        can_lower_scalar_h1_stiffness_localform_v<ExprT>,
        "lower_scalar_h1_stiffness_localform only supports scalar coefficient * inner(gradt(u),grad(v)) local forms" );
    using lowered_expr_type = LoweredScalarStiffnessExpr<
        typename factors::test_space_type,
        typename factors::trial_space_type,
        typename factors::coefficient_expr_type>;
    return vf::expr( lowered_expr_type(
        factors::coefficient_expr( expr ), expr.polynomialOrder() ) );
}

/**
 * @brief Lower `gradt(u)*K*trans(grad(v))` to the shared local tensor shell.
 * @param expr canonical anisotropic scalar H1 expression
 * @return lowered expression carrying the original tensor coefficient
 */
template<typename ExprT>
auto lower_anisotropic_scalar_h1_stiffness_localform( ExprT const& expr )
{
    using factors = anisotropic_scalar_h1_stiffness_factors<ExprT>;
    static_assert(
        can_lower_anisotropic_scalar_h1_stiffness_localform_v<ExprT>,
        "lower_anisotropic_scalar_h1_stiffness_localform expects gradt(u)*K*trans(grad(v)) with a square tensor coefficient" );
    using lowered_expr_type = LoweredScalarBilinearExpr<
        typename factors::test_space_type,
        typename factors::trial_space_type,
        typename factors::coefficient_expr_type,
        anisotropic_scalar_h1_stiffness_contraction>;
    return vf::expr( lowered_expr_type(
        factors::coefficient_expr( expr ), expr.polynomialOrder() ) );
}

/**
 * @brief Lower a canonical operator pair with a selected contraction policy.
 *
 * @tparam ContractionType CPU/Eigen mathematical contraction strategy
 * @tparam ExprT recognized canonical bilinear expression
 * @param expr original DSEL expression
 * @return lowered expression using the shared local-tensor shell
 */
template<typename ContractionType, typename ExprT>
auto lower_canonical_bilinear_localform( ExprT const& expr )
{
    using factors = canonical_bilinear_factors<ExprT>;
    static_assert( factors::supported && factors::test_count == 1 &&
                   factors::trial_count == 1,
                   "lower_canonical_bilinear_localform requires one canonical test and trial operator" );
    using lowered_expr_type = LoweredScalarBilinearExpr<
        typename factors::test_space_type,
        typename factors::trial_space_type,
        typename factors::coefficient_expr_type,
        ContractionType>;
    return vf::expr( lowered_expr_type(
        factors::coefficient_expr( expr ), expr.polynomialOrder() ) );
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
