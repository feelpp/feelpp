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
