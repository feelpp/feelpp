/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file localform_base.hpp
   \brief Shared local-form lowering utilities and tags.
 */
#ifndef FEELPP_VF_DETAIL_LOCALFORM_BASE_HPP
#define FEELPP_VF_DETAIL_LOCALFORM_BASE_HPP 1

#include <type_traits>
#include <utility>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/operators.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

template<typename T>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

template<typename T1, typename T2>
using first_non_void_t = std::conditional_t<!std::is_void_v<T1>, T1, T2>;

template<typename T1, typename T2>
inline constexpr bool same_or_void_v =
    std::is_void_v<T1> || std::is_void_v<T2> || std::is_same_v<T1, T2>;

template<typename TestSpaceType, typename TrialSpaceType, typename CoefficientExprType>
class LoweredScalarMassExpr;

template<typename SpaceType, typename CoefficientExprType>
class LoweredScalarSourceExpr;

template<typename T>
struct is_lowered_localform_expr : std::false_type {};

template<typename TestSpaceType, typename TrialSpaceType, typename CoefficientExprType>
struct is_lowered_localform_expr<LoweredScalarMassExpr<TestSpaceType, TrialSpaceType, CoefficientExprType>> : std::true_type {};

template<typename SpaceType, typename CoefficientExprType>
struct is_lowered_localform_expr<LoweredScalarSourceExpr<SpaceType, CoefficientExprType>> : std::true_type {};

template<typename T>
struct is_lowered_localform : std::false_type {};

template<typename ExprT>
struct is_lowered_localform<Expr<ExprT>> : is_lowered_localform_expr<ExprT> {};

template<typename T>
inline constexpr bool is_lowered_localform_v = is_lowered_localform<remove_cvref_t<T>>::value;

enum class localform_scalar_expr_kind
{
    unsupported,
    scalar_leaf,
    scalar_value_leaf,
    geometry_scalar_leaf,
    unary_composite,
    binary_composite
};

template<typename ExprT>
struct is_localform_cst_expr : std::false_type {};

template<typename T>
struct is_localform_cst_expr<Expr<Cst<T>>> : std::true_type {};

template<typename ExprT>
struct is_localform_value_id_expr : std::false_type {};

template<typename Element>
struct is_localform_value_id_expr<Expr<OpId<Element, __VALUE>>> : std::true_type {};

template<typename T, typename = void>
struct has_localform_unary_subexpression : std::false_type {};

template<typename T>
struct has_localform_unary_subexpression<T,
                                         std::void_t<decltype( std::declval<remove_cvref_t<T> const&>().expression().expression() )>>
    : std::true_type {};

template<typename T, typename = void>
struct has_localform_binary_subexpressions : std::false_type {};

template<typename T>
struct has_localform_binary_subexpressions<T,
                                           std::void_t<decltype( std::declval<remove_cvref_t<T> const&>().expression().expression1() ),
                                                       decltype( std::declval<remove_cvref_t<T> const&>().expression().expression2() )>>
    : std::true_type {};

template<typename T, typename = void>
struct has_localform_left_right_subexpressions : std::false_type {};

template<typename T>
struct has_localform_left_right_subexpressions<T,
                                               std::void_t<decltype( std::declval<remove_cvref_t<T> const&>().expression().left() ),
                                                           decltype( std::declval<remove_cvref_t<T> const&>().expression().right() )>>
    : std::true_type {};

template<typename T, typename = void>
struct is_localform_scalar_evaluate_type : std::false_type {};

template<typename T, typename = void>
struct has_localform_fixed_matrix_extents : std::false_type {};

template<typename T>
struct has_localform_fixed_matrix_extents<T, std::void_t<decltype( remove_cvref_t<T>::RowsAtCompileTime ),
                                                          decltype( remove_cvref_t<T>::ColsAtCompileTime )>>
    : std::true_type {};

template<typename T>
struct is_localform_scalar_evaluate_type<T, std::void_t<typename remove_cvref_t<T>::evaluate_type>>
{
    using evaluate_type = typename remove_cvref_t<T>::evaluate_type;
    static inline const bool value =
        std::is_arithmetic_v<evaluate_type> ||
        ( has_localform_fixed_matrix_extents<evaluate_type>::value &&
          evaluate_type::RowsAtCompileTime == 1 &&
          evaluate_type::ColsAtCompileTime == 1 );
};

} // namespace detail
} // namespace vf
} // namespace Feel

#endif
