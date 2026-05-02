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
   \file localform_match_impl.hpp
   \brief Matching and factorization for lowered scalar local forms.
 */
#ifndef FEELPP_VF_DETAIL_LOCALFORM_MATCH_IMPL_HPP
#define FEELPP_VF_DETAIL_LOCALFORM_MATCH_IMPL_HPP 1

#include <feel/feelvf/detail/localform_base.hpp>
#include <feel/feelvf/detail/localform_coeff_ir.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/operations.hpp>

namespace Feel
{
namespace vf
{
template<typename ExprT1>
class __Chi__;
}
}

namespace Feel
{
namespace vf
{
namespace detail
{

template<typename ExprT>
struct is_localform_geometry_scalar_expr : std::false_type {};

#define FEELPP_LOCALFORM_GEOMETRY_SCALAR_LEAFS(_) \
    _( GDNx ) \
    _( GDNy ) \
    _( GDNz ) \
    _( GDnormalNorm ) \
    _( GDTx ) \
    _( GDTy ) \
    _( GDTz ) \
    _( GDDetJ ) \
    _( GDPx ) \
    _( GDPy ) \
    _( GDPz ) \
    _( GDCx ) \
    _( GDCy ) \
    _( GDCz ) \
    _( GDH ) \
    _( GDHMin ) \
    _( GDHFace ) \
    _( GDMeas ) \
    _( GDMeasPEN ) \
    _( GDNPEN ) \
    _( GDHMeasFace ) \
    _( GDEid ) \
    _( GDEmarker ) \
    _( GDFmarker ) \
    _( GDEmarker2 ) \
    _( GDEPid )

#define FEELPP_LOCALFORM_DECLARE_GEOMETRY_SCALAR_EXPR(ExprName) \
    template<> \
    struct is_localform_geometry_scalar_expr<Expr<ExprName>> : std::true_type {};

FEELPP_LOCALFORM_GEOMETRY_SCALAR_LEAFS( FEELPP_LOCALFORM_DECLARE_GEOMETRY_SCALAR_EXPR )

#undef FEELPP_LOCALFORM_DECLARE_GEOMETRY_SCALAR_EXPR

template<bool IsSupported, typename SupportedType>
using localform_supported_coeff_expr_t = std::conditional_t<IsSupported, SupportedType, void>;

template<typename ExprT, typename Enable = void>
struct scalar_localform_coefficient_expr
{
    using expr_type = remove_cvref_t<ExprT>;
    using coefficient_expr_type = void;
    static constexpr auto kind = localform_scalar_expr_kind::unsupported;
    static inline const bool supported = false;

    static coefficient_expr_type coefficient_expr( expr_type const& ) = delete;
};

template<typename T>
struct scalar_localform_coefficient_expr<Expr<Cst<T>>, void>
{
    using expr_type = Expr<Cst<T>>;
    using coefficient_expr_type = LocalformCoeffConstant<expr_type>;
    static constexpr auto kind = localform_scalar_expr_kind::scalar_leaf;
    static inline const bool supported = true;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( expr );
    }
};

template<typename Element>
struct scalar_localform_coefficient_expr<Expr<OpId<Element, __VALUE>>, void>
{
    using expr_type = Expr<OpId<Element, __VALUE>>;
    using space_type = typename remove_cvref_t<Element>::functionspace_type;
    using coefficient_expr_type = LocalformCoeffScalarValue<expr_type>;
    static constexpr auto kind = localform_scalar_expr_kind::scalar_value_leaf;
    static inline const bool supported = space_type::nSpaces == 1 && space_type::is_scalar;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( expr );
    }
};

template<typename ExprT>
struct scalar_localform_coefficient_expr<ExprT,
                                         std::enable_if_t<is_localform_geometry_scalar_expr<remove_cvref_t<ExprT>>::value>>
{
    using expr_type = remove_cvref_t<ExprT>;
    using coefficient_expr_type = LocalformCoeffGeometry<expr_type>;
    static constexpr auto kind = localform_scalar_expr_kind::geometry_scalar_leaf;
    static inline const bool supported = true;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( expr );
    }
};

template<typename InnerExpr>
struct scalar_localform_coefficient_expr<Expr<__Chi__<InnerExpr>>, void>
{
    using expr_type = Expr<__Chi__<InnerExpr>>;
    using inner_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().expression() )>;
    using inner_support = scalar_localform_coefficient_expr<inner_expr_type>;
    using coefficient_expr_type =
        localform_supported_coeff_expr_t<inner_support::supported,
                                         LocalformCoeffChi<typename inner_support::coefficient_expr_type>>;
    static constexpr auto kind = localform_scalar_expr_kind::unary_composite;
    static inline const bool supported = inner_support::supported;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( inner_support::coefficient_expr( expr.expression().expression() ) );
    }
};

template<typename L, typename R>
struct scalar_localform_coefficient_expr<Expr<vf_add<L, R>>, void>
{
    using expr_type = Expr<vf_add<L, R>>;
    using left_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().right() )>;
    using left_support = scalar_localform_coefficient_expr<left_expr_type>;
    using right_support = scalar_localform_coefficient_expr<right_expr_type>;
    using coefficient_expr_type =
        localform_supported_coeff_expr_t<left_support::supported && right_support::supported,
                                         LocalformCoeffAdd<typename left_support::coefficient_expr_type,
                                                           typename right_support::coefficient_expr_type>>;
    static constexpr auto kind = localform_scalar_expr_kind::binary_composite;
    static inline const bool supported = left_support::supported && right_support::supported;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( left_support::coefficient_expr( expr.expression().left() ),
                                      right_support::coefficient_expr( expr.expression().right() ) );
    }
};

template<typename L, typename R>
struct scalar_localform_coefficient_expr<Expr<vf_mul<L, R>>, void>
{
    using expr_type = Expr<vf_mul<L, R>>;
    using left_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().right() )>;
    using left_support = scalar_localform_coefficient_expr<left_expr_type>;
    using right_support = scalar_localform_coefficient_expr<right_expr_type>;
    using coefficient_expr_type =
        localform_supported_coeff_expr_t<left_support::supported && right_support::supported,
                                         LocalformCoeffMul<typename left_support::coefficient_expr_type,
                                                           typename right_support::coefficient_expr_type>>;
    static constexpr auto kind = localform_scalar_expr_kind::binary_composite;
    static inline const bool supported = left_support::supported && right_support::supported;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( left_support::coefficient_expr( expr.expression().left() ),
                                      right_support::coefficient_expr( expr.expression().right() ) );
    }
};

template<typename L, typename R>
struct scalar_localform_coefficient_expr<Expr<vf_greater<L, R>>, void>
{
    using expr_type = Expr<vf_greater<L, R>>;
    using left_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().right() )>;
    using left_support = scalar_localform_coefficient_expr<left_expr_type>;
    using right_support = scalar_localform_coefficient_expr<right_expr_type>;
    using coefficient_expr_type =
        localform_supported_coeff_expr_t<left_support::supported && right_support::supported,
                                         LocalformCoeffGreater<typename left_support::coefficient_expr_type,
                                                               typename right_support::coefficient_expr_type>>;
    static constexpr auto kind = localform_scalar_expr_kind::binary_composite;
    static inline const bool supported = left_support::supported && right_support::supported;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( left_support::coefficient_expr( expr.expression().left() ),
                                      right_support::coefficient_expr( expr.expression().right() ) );
    }
};

template<typename ExprT>
struct scalar_localform_factors
{
    using expr_type = remove_cvref_t<ExprT>;
    using coefficient_support = scalar_localform_coefficient_expr<expr_type>;
    static inline const bool supported = coefficient_support::supported;
    static inline const int test_count = 0;
    static inline const int trial_count = 0;
    using coefficient_expr_type = typename coefficient_support::coefficient_expr_type;
    using test_space_type = void;
    using trial_space_type = void;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_support::coefficient_expr( expr );
    }
};

template<typename Element>
struct scalar_localform_factors<Expr<OpId<Element, __TEST>>>
{
    using expr_type = Expr<OpId<Element, __TEST>>;
    using space_type = typename remove_cvref_t<Element>::functionspace_type;
    using coefficient_expr_type = LocalformCoeffConstant<Expr<Cst<typename expr_type::value_type>>>;
    static inline const bool supported = space_type::nSpaces == 1 && space_type::is_scalar;
    static inline const int test_count = 1;
    static inline const int trial_count = 0;
    using test_space_type = space_type;
    using trial_space_type = void;

    static coefficient_expr_type coefficient_expr( expr_type const& )
    {
        return coefficient_expr_type( vf::cst( typename expr_type::value_type( 1 ) ) );
    }
};

template<typename Element>
struct scalar_localform_factors<Expr<OpId<Element, __TRIAL>>>
{
    using expr_type = Expr<OpId<Element, __TRIAL>>;
    using space_type = typename remove_cvref_t<Element>::functionspace_type;
    using coefficient_expr_type = LocalformCoeffConstant<Expr<Cst<typename expr_type::value_type>>>;
    static inline const bool supported = space_type::nSpaces == 1 && space_type::is_scalar;
    static inline const int test_count = 0;
    static inline const int trial_count = 1;
    using test_space_type = void;
    using trial_space_type = space_type;

    static coefficient_expr_type coefficient_expr( expr_type const& )
    {
        return coefficient_expr_type( vf::cst( typename expr_type::value_type( 1 ) ) );
    }
};

template<typename L, typename R>
struct scalar_localform_factors<Expr<vf_mul<L, R>>>
{
    using expr_type = Expr<vf_mul<L, R>>;
    using left_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype( std::declval<expr_type const&>().expression().right() )>;
    using left_factors = scalar_localform_factors<left_expr_type>;
    using right_factors = scalar_localform_factors<right_expr_type>;
    using coefficient_expr_type =
        localform_supported_coeff_expr_t<left_factors::supported &&
                                         right_factors::supported &&
                                         same_or_void_v<typename left_factors::test_space_type, typename right_factors::test_space_type> &&
                                         same_or_void_v<typename left_factors::trial_space_type, typename right_factors::trial_space_type> &&
                                         ( left_factors::test_count + right_factors::test_count <= 1 ) &&
                                         ( left_factors::trial_count + right_factors::trial_count <= 1 ),
                                         LocalformCoeffMul<typename left_factors::coefficient_expr_type,
                                                           typename right_factors::coefficient_expr_type>>;
    static inline const bool supported =
        left_factors::supported &&
        right_factors::supported &&
        same_or_void_v<typename left_factors::test_space_type, typename right_factors::test_space_type> &&
        same_or_void_v<typename left_factors::trial_space_type, typename right_factors::trial_space_type> &&
        ( left_factors::test_count + right_factors::test_count <= 1 ) &&
        ( left_factors::trial_count + right_factors::trial_count <= 1 );
    static inline const int test_count = left_factors::test_count + right_factors::test_count;
    static inline const int trial_count = left_factors::trial_count + right_factors::trial_count;
    using test_space_type = first_non_void_t<typename left_factors::test_space_type,
                                             typename right_factors::test_space_type>;
    using trial_space_type = first_non_void_t<typename left_factors::trial_space_type,
                                              typename right_factors::trial_space_type>;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type( left_factors::coefficient_expr( expr.expression().left() ),
                                      right_factors::coefficient_expr( expr.expression().right() ) );
    }
};

template<typename ExprT>
inline constexpr bool can_lower_scalar_bilinear_localform_v =
    !is_lowered_localform_v<ExprT> &&
    scalar_localform_factors<ExprT>::supported &&
    scalar_localform_factors<ExprT>::test_count == 1 &&
    scalar_localform_factors<ExprT>::trial_count == 1;

template<typename ExprT>
inline constexpr bool can_lower_scalar_linear_localform_v =
    !is_lowered_localform_v<ExprT> &&
    scalar_localform_factors<ExprT>::supported &&
    scalar_localform_factors<ExprT>::test_count == 1 &&
    scalar_localform_factors<ExprT>::trial_count == 0;

template<typename ExprT, typename FormType>
inline constexpr bool can_lower_scalar_bilinear_localform_for_form_v =
    can_lower_scalar_bilinear_localform_v<ExprT> &&
    FormType::test_space_type::nSpaces == 1 &&
    FormType::trial_space_type::nSpaces == 1 &&
    FormType::test_space_type::is_scalar &&
    FormType::trial_space_type::is_scalar &&
    std::is_same_v<typename scalar_localform_factors<ExprT>::test_space_type,
                   typename FormType::test_space_type> &&
    std::is_same_v<typename scalar_localform_factors<ExprT>::trial_space_type,
                   typename FormType::trial_space_type>;

template<typename ExprT, typename FormType>
inline constexpr bool can_lower_scalar_linear_localform_for_form_v =
    can_lower_scalar_linear_localform_v<ExprT> &&
    FormType::test_space_type::nSpaces == 1 &&
    FormType::test_space_type::is_scalar &&
    std::is_same_v<typename scalar_localform_factors<ExprT>::test_space_type,
                   typename FormType::test_space_type>;

} // namespace detail
} // namespace vf
} // namespace Feel

#undef FEELPP_LOCALFORM_GEOMETRY_SCALAR_LEAFS

#endif
