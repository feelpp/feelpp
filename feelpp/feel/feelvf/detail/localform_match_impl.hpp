/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
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
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/trans.hpp>

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

/** @brief Unsupported scalar H1 gradient factor. */
template<typename ExprT>
struct scalar_h1_gradient_factor
{
    static inline const bool supported = false;
    static inline const int test_count = 0;
    static inline const int trial_count = 0;
    using test_space_type = void;
    using trial_space_type = void;
};

/** @brief Scalar test-function gradient factor `grad(v)`. */
template<typename Element>
struct scalar_h1_gradient_factor<Expr<OpGrad<Element, __TEST>>>
{
    using space_type = typename remove_cvref_t<Element>::functionspace_type;
    static inline const bool supported = space_type::nSpaces == 1 &&
                                         space_type::is_scalar;
    static inline const int test_count = 1;
    static inline const int trial_count = 0;
    using test_space_type = space_type;
    using trial_space_type = void;
};

/** @brief Scalar trial-function gradient factor `gradt(u)`. */
template<typename Element>
struct scalar_h1_gradient_factor<Expr<OpGrad<Element, __TRIAL>>>
{
    using space_type = typename remove_cvref_t<Element>::functionspace_type;
    static inline const bool supported = space_type::nSpaces == 1 &&
                                         space_type::is_scalar;
    static inline const int test_count = 0;
    static inline const int trial_count = 1;
    using test_space_type = void;
    using trial_space_type = space_type;
};

/**
 * @brief Factors of a scalar isotropic H1 stiffness expression.
 *
 * The recognized mathematical form is a supported scalar coefficient times
 * `inner(gradt(u),grad(v))`; multiplication may be associated arbitrarily.
 */
template<typename ExprT, typename Enable = void>
struct scalar_h1_stiffness_factors
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

/** @brief Factorization of `inner(gradt(u),grad(v))`. */
template<typename L, typename R>
struct scalar_h1_stiffness_factors<Expr<Product<L, R, 1, NONE>>, void>
{
    using expr_type = Expr<Product<L, R, 1, NONE>>;
    using left_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().right() )>;
    using left_gradient = scalar_h1_gradient_factor<left_expr_type>;
    using right_gradient = scalar_h1_gradient_factor<right_expr_type>;
    using coefficient_expr_type =
        LocalformCoeffConstant<Expr<Cst<typename expr_type::value_type>>>;
    static inline const bool supported =
        left_gradient::supported && right_gradient::supported &&
        same_or_void_v<typename left_gradient::test_space_type,
                       typename right_gradient::test_space_type> &&
        same_or_void_v<typename left_gradient::trial_space_type,
                       typename right_gradient::trial_space_type> &&
        left_gradient::test_count + right_gradient::test_count == 1 &&
        left_gradient::trial_count + right_gradient::trial_count == 1;
    static inline const int test_count = left_gradient::test_count +
                                         right_gradient::test_count;
    static inline const int trial_count = left_gradient::trial_count +
                                          right_gradient::trial_count;
    using test_space_type = first_non_void_t<
        typename left_gradient::test_space_type,
        typename right_gradient::test_space_type>;
    using trial_space_type = first_non_void_t<
        typename left_gradient::trial_space_type,
        typename right_gradient::trial_space_type>;

    static coefficient_expr_type coefficient_expr( expr_type const& )
    {
        return coefficient_expr_type(
            vf::cst( typename expr_type::value_type( 1 ) ) );
    }
};

/** @brief Factorization of a multiplied scalar H1 stiffness expression. */
template<typename L, typename R>
struct scalar_h1_stiffness_factors<Expr<vf_mul<L, R>>, void>
{
    using expr_type = Expr<vf_mul<L, R>>;
    using left_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().right() )>;
    using left_factors = scalar_h1_stiffness_factors<left_expr_type>;
    using right_factors = scalar_h1_stiffness_factors<right_expr_type>;
    using coefficient_expr_type = localform_supported_coeff_expr_t<
        left_factors::supported && right_factors::supported &&
        same_or_void_v<typename left_factors::test_space_type,
                       typename right_factors::test_space_type> &&
        same_or_void_v<typename left_factors::trial_space_type,
                       typename right_factors::trial_space_type> &&
        left_factors::test_count + right_factors::test_count <= 1 &&
        left_factors::trial_count + right_factors::trial_count <= 1,
        LocalformCoeffMul<typename left_factors::coefficient_expr_type,
                          typename right_factors::coefficient_expr_type>>;
    static inline const bool supported =
        left_factors::supported && right_factors::supported &&
        same_or_void_v<typename left_factors::test_space_type,
                       typename right_factors::test_space_type> &&
        same_or_void_v<typename left_factors::trial_space_type,
                       typename right_factors::trial_space_type> &&
        left_factors::test_count + right_factors::test_count <= 1 &&
        left_factors::trial_count + right_factors::trial_count <= 1;
    static inline const int test_count = left_factors::test_count +
                                         right_factors::test_count;
    static inline const int trial_count = left_factors::trial_count +
                                          right_factors::trial_count;
    using test_space_type = first_non_void_t<
        typename left_factors::test_space_type,
        typename right_factors::test_space_type>;
    using trial_space_type = first_non_void_t<
        typename left_factors::trial_space_type,
        typename right_factors::trial_space_type>;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type(
            left_factors::coefficient_expr( expr.expression().left() ),
            right_factors::coefficient_expr( expr.expression().right() ) );
    }
};

/** @brief Canonical FE operator carried by one test or trial factor. */
enum class localform_basis_operator
{
    none,
    value,
    gradient,
    divergence,
    curl,
    unsupported
};

/** @brief Compile-time role and operator of one canonical FE factor. */
template<typename ExprT>
struct canonical_basis_factor
{
    static inline const bool supported = false;
    static inline const int test_count = 0;
    static inline const int trial_count = 0;
    static constexpr auto test_operator = localform_basis_operator::none;
    static constexpr auto trial_operator = localform_basis_operator::none;
    using test_space_type = void;
    using trial_space_type = void;
};

#define FEELPP_CANONICAL_BASIS_FACTOR(OP, ROLE, TEST_COUNT, TRIAL_COUNT, TEST_OP, TRIAL_OP, TEST_SPACE, TRIAL_SPACE) \
    template<typename Element> \
    struct canonical_basis_factor<Expr<OP<Element, ROLE>>> \
    { \
        using expr_type = Expr<OP<Element, ROLE>>; \
        using space_type = typename remove_cvref_t<Element>::functionspace_type; \
        using coefficient_expr_type = LocalformCoeffConstant<Expr<Cst<typename expr_type::value_type>>>; \
        static inline const bool supported = space_type::nSpaces == 1; \
        static inline const int test_count = TEST_COUNT; \
        static inline const int trial_count = TRIAL_COUNT; \
        static constexpr auto test_operator = localform_basis_operator::TEST_OP; \
        static constexpr auto trial_operator = localform_basis_operator::TRIAL_OP; \
        using test_space_type = TEST_SPACE; \
        using trial_space_type = TRIAL_SPACE; \
        static coefficient_expr_type coefficient_expr( expr_type const& ) \
        { \
            return coefficient_expr_type( vf::cst( typename expr_type::value_type( 1 ) ) ); \
        } \
    };

FEELPP_CANONICAL_BASIS_FACTOR( OpId, __TEST, 1, 0, value, none,
                               space_type, void )
FEELPP_CANONICAL_BASIS_FACTOR( OpId, __TRIAL, 0, 1, none, value,
                               void, space_type )
FEELPP_CANONICAL_BASIS_FACTOR( OpGrad, __TEST, 1, 0, gradient, none,
                               space_type, void )
FEELPP_CANONICAL_BASIS_FACTOR( OpGrad, __TRIAL, 0, 1, none, gradient,
                               void, space_type )
FEELPP_CANONICAL_BASIS_FACTOR( OpDiv, __TEST, 1, 0, divergence, none,
                               space_type, void )
FEELPP_CANONICAL_BASIS_FACTOR( OpDiv, __TRIAL, 0, 1, none, divergence,
                               void, space_type )
FEELPP_CANONICAL_BASIS_FACTOR( OpCurl, __TEST, 1, 0, curl, none,
                               space_type, void )
FEELPP_CANONICAL_BASIS_FACTOR( OpCurl, __TRIAL, 0, 1, none, curl,
                               void, space_type )

#undef FEELPP_CANONICAL_BASIS_FACTOR

/** @return the only operator in two factors, or unsupported on collision. */
consteval localform_basis_operator
mergeBasisOperator( localform_basis_operator left,
                    localform_basis_operator right )
{
    if ( left == localform_basis_operator::none )
        return right;
    if ( right == localform_basis_operator::none )
        return left;
    return localform_basis_operator::unsupported;
}

/**
 * @brief Shared factorization for canonical value/differential bilinear forms.
 */
template<typename ExprT, typename Enable = void>
struct canonical_bilinear_factors
{
    using expr_type = remove_cvref_t<ExprT>;
    using coefficient_support = scalar_localform_coefficient_expr<expr_type>;
    static inline const bool supported = coefficient_support::supported;
    static inline const int test_count = 0;
    static inline const int trial_count = 0;
    static constexpr auto test_operator = localform_basis_operator::none;
    static constexpr auto trial_operator = localform_basis_operator::none;
    using coefficient_expr_type = typename coefficient_support::coefficient_expr_type;
    using test_space_type = void;
    using trial_space_type = void;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_support::coefficient_expr( expr );
    }
};

/** @brief Factorization of an inner product between two FE operators. */
template<typename L, typename R>
struct canonical_bilinear_factors<Expr<Product<L, R, 1, NONE>>, void>
{
    using expr_type = Expr<Product<L, R, 1, NONE>>;
    using left_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().right() )>;
    using left_factor = canonical_basis_factor<left_expr_type>;
    using right_factor = canonical_basis_factor<right_expr_type>;
    static constexpr auto test_operator = mergeBasisOperator(
        left_factor::test_operator, right_factor::test_operator );
    static constexpr auto trial_operator = mergeBasisOperator(
        left_factor::trial_operator, right_factor::trial_operator );
    using coefficient_expr_type =
        LocalformCoeffConstant<Expr<Cst<typename expr_type::value_type>>>;
    static inline const bool supported =
        left_factor::supported && right_factor::supported &&
        test_operator != localform_basis_operator::unsupported &&
        trial_operator != localform_basis_operator::unsupported &&
        left_factor::test_count + right_factor::test_count == 1 &&
        left_factor::trial_count + right_factor::trial_count == 1;
    static inline const int test_count = left_factor::test_count +
                                         right_factor::test_count;
    static inline const int trial_count = left_factor::trial_count +
                                          right_factor::trial_count;
    using test_space_type = first_non_void_t<
        typename left_factor::test_space_type,
        typename right_factor::test_space_type>;
    using trial_space_type = first_non_void_t<
        typename left_factor::trial_space_type,
        typename right_factor::trial_space_type>;

    static coefficient_expr_type coefficient_expr( expr_type const& )
    {
        return coefficient_expr_type(
            vf::cst( typename expr_type::value_type( 1 ) ) );
    }
};

/** @brief Factorization of multiplied coefficients and canonical FE factors. */
template<typename L, typename R>
struct canonical_bilinear_factors<Expr<vf_mul<L, R>>, void>
{
    using expr_type = Expr<vf_mul<L, R>>;
    using left_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().right() )>;
    using left_factors = std::conditional_t<
        canonical_basis_factor<left_expr_type>::supported,
        canonical_basis_factor<left_expr_type>,
        canonical_bilinear_factors<left_expr_type>>;
    using right_factors = std::conditional_t<
        canonical_basis_factor<right_expr_type>::supported,
        canonical_basis_factor<right_expr_type>,
        canonical_bilinear_factors<right_expr_type>>;
    static constexpr auto test_operator = mergeBasisOperator(
        left_factors::test_operator, right_factors::test_operator );
    static constexpr auto trial_operator = mergeBasisOperator(
        left_factors::trial_operator, right_factors::trial_operator );
    static inline const bool structurally_supported =
        left_factors::supported && right_factors::supported &&
        test_operator != localform_basis_operator::unsupported &&
        trial_operator != localform_basis_operator::unsupported &&
        left_factors::test_count + right_factors::test_count <= 1 &&
        left_factors::trial_count + right_factors::trial_count <= 1 &&
        same_or_void_v<typename left_factors::test_space_type,
                       typename right_factors::test_space_type> &&
        same_or_void_v<typename left_factors::trial_space_type,
                       typename right_factors::trial_space_type>;
    using coefficient_expr_type = localform_supported_coeff_expr_t<
        structurally_supported,
        LocalformCoeffMul<typename left_factors::coefficient_expr_type,
                          typename right_factors::coefficient_expr_type>>;
    static inline const bool supported = structurally_supported;
    static inline const int test_count = left_factors::test_count +
                                         right_factors::test_count;
    static inline const int trial_count = left_factors::trial_count +
                                          right_factors::trial_count;
    using test_space_type = first_non_void_t<
        typename left_factors::test_space_type,
        typename right_factors::test_space_type>;
    using trial_space_type = first_non_void_t<
        typename left_factors::trial_space_type,
        typename right_factors::trial_space_type>;

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return coefficient_expr_type(
            left_factors::coefficient_expr( expr.expression().left() ),
            right_factors::coefficient_expr( expr.expression().right() ) );
    }
};

/** @brief Canonical operator-pair predicate. */
template<typename ExprT, localform_basis_operator TestOperator,
         localform_basis_operator TrialOperator>
inline constexpr bool can_lower_canonical_bilinear_localform_v =
    !is_lowered_localform_v<ExprT> &&
    canonical_bilinear_factors<ExprT>::supported &&
    canonical_bilinear_factors<ExprT>::test_count == 1 &&
    canonical_bilinear_factors<ExprT>::trial_count == 1 &&
    canonical_bilinear_factors<ExprT>::test_operator == TestOperator &&
    canonical_bilinear_factors<ExprT>::trial_operator == TrialOperator;

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

/** @brief True when an expression is a canonical scalar H1 stiffness form. */
template<typename ExprT>
inline constexpr bool can_lower_scalar_h1_stiffness_localform_v =
    !is_lowered_localform_v<ExprT> &&
    scalar_h1_stiffness_factors<ExprT>::supported &&
    scalar_h1_stiffness_factors<ExprT>::test_count == 1 &&
    scalar_h1_stiffness_factors<ExprT>::trial_count == 1;

/** @brief Unsupported trial-gradient/tensor product. */
template<typename ExprT>
struct anisotropic_h1_trial_tensor_factor
{
    static constexpr bool supported = false;
    using coefficient_expr_type = void;
    using trial_space_type = void;
};

/** @brief Factorization of `gradt(u)*K` for a matrix coefficient `K`. */
template<typename L, typename R>
struct anisotropic_h1_trial_tensor_factor<Expr<vf_mul<L, R>>>
{
    using expr_type = Expr<vf_mul<L, R>>;
    using gradient_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().left() )>;
    using coefficient_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().right() )>;
    using gradient_factor = scalar_h1_gradient_factor<gradient_expr_type>;
    using trial_space_type = typename gradient_factor::trial_space_type;
    static constexpr bool supported = []()
    {
        if constexpr (
            gradient_factor::supported && gradient_factor::test_count == 0 &&
            gradient_factor::trial_count == 1 &&
            !std::is_same_v<trial_space_type, void> &&
            requires {
                typename coefficient_expr_type::evaluate_type;
                typename coefficient_expr_type::test_basis;
                typename coefficient_expr_type::trial_basis;
                coefficient_expr_type::evaluate_type::RowsAtCompileTime;
                coefficient_expr_type::evaluate_type::ColsAtCompileTime;
            } )
        {
            return coefficient_expr_type::evaluate_type::RowsAtCompileTime ==
                       trial_space_type::nDim &&
                   coefficient_expr_type::evaluate_type::ColsAtCompileTime ==
                       trial_space_type::nDim &&
                   std::is_same_v<typename coefficient_expr_type::test_basis,
                                  std::nullptr_t> &&
                   std::is_same_v<typename coefficient_expr_type::trial_basis,
                                  std::nullptr_t>;
        }
        else
            return false;
    }();

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return expr.expression().right();
    }
};

/** @brief Unsupported transposed scalar test gradient. */
template<typename ExprT>
struct anisotropic_h1_test_gradient_factor
{
    static constexpr bool supported = false;
    using test_space_type = void;
};

/** @brief Factorization of `trans(grad(v))`. */
template<typename InnerExpr>
struct anisotropic_h1_test_gradient_factor<Expr<Trans<InnerExpr>>>
{
    using gradient_expr_type = remove_cvref_t<InnerExpr>;
    using gradient_factor = scalar_h1_gradient_factor<gradient_expr_type>;
    using test_space_type = typename gradient_factor::test_space_type;
    static constexpr bool supported =
        gradient_factor::supported && gradient_factor::test_count == 1 &&
        gradient_factor::trial_count == 0;
};

/** @brief Unsupported anisotropic scalar H1 expression. */
template<typename ExprT>
struct anisotropic_scalar_h1_stiffness_factors
{
    static constexpr bool supported = false;
    using coefficient_expr_type = void;
    using test_space_type = void;
    using trial_space_type = void;
};

/**
 * @brief Factorization of the canonical `gradt(u)*K*trans(grad(v))` form.
 */
template<typename L, typename R>
struct anisotropic_scalar_h1_stiffness_factors<Expr<vf_mul<L, R>>>
{
    using expr_type = Expr<vf_mul<L, R>>;
    using left_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().left() )>;
    using right_expr_type = remove_cvref_t<decltype(
        std::declval<expr_type const&>().expression().right() )>;
    using trial_tensor = anisotropic_h1_trial_tensor_factor<left_expr_type>;
    using test_gradient = anisotropic_h1_test_gradient_factor<right_expr_type>;
    using coefficient_expr_type = typename trial_tensor::coefficient_expr_type;
    using test_space_type = typename test_gradient::test_space_type;
    using trial_space_type = typename trial_tensor::trial_space_type;
    static constexpr bool supported = []()
    {
        if constexpr (
            trial_tensor::supported && test_gradient::supported &&
            !std::is_same_v<test_space_type, void> &&
            !std::is_same_v<trial_space_type, void> )
            return test_space_type::nDim == trial_space_type::nDim;
        else
            return false;
    }();

    static coefficient_expr_type coefficient_expr( expr_type const& expr )
    {
        return trial_tensor::coefficient_expr( expr.expression().left() );
    }
};

/** @brief True for canonical tensor-diffusion syntax. */
template<typename ExprT>
inline constexpr bool can_lower_anisotropic_scalar_h1_stiffness_localform_v =
    !is_lowered_localform_v<ExprT> &&
    anisotropic_scalar_h1_stiffness_factors<ExprT>::supported;

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

/** @brief True when a stiffness expression and bilinear form spaces agree. */
template<typename ExprT, typename FormType>
inline constexpr bool can_lower_scalar_h1_stiffness_localform_for_form_v =
    can_lower_scalar_h1_stiffness_localform_v<ExprT> &&
    FormType::test_space_type::nSpaces == 1 &&
    FormType::trial_space_type::nSpaces == 1 &&
    FormType::test_space_type::is_scalar &&
    FormType::trial_space_type::is_scalar &&
    std::is_same_v<
        typename scalar_h1_stiffness_factors<ExprT>::test_space_type,
        typename FormType::test_space_type> &&
    std::is_same_v<
        typename scalar_h1_stiffness_factors<ExprT>::trial_space_type,
        typename FormType::trial_space_type>;

/** @brief True when anisotropic H1 syntax and form spaces agree. */
template<typename ExprT, typename FormType>
inline constexpr bool
can_lower_anisotropic_scalar_h1_stiffness_localform_for_form_v =
    can_lower_anisotropic_scalar_h1_stiffness_localform_v<ExprT> &&
    FormType::test_space_type::nSpaces == 1 &&
    FormType::trial_space_type::nSpaces == 1 &&
    FormType::test_space_type::is_scalar &&
    FormType::trial_space_type::is_scalar &&
    std::is_same_v<
        typename anisotropic_scalar_h1_stiffness_factors<
            ExprT>::test_space_type,
        typename FormType::test_space_type> &&
    std::is_same_v<
        typename anisotropic_scalar_h1_stiffness_factors<
            ExprT>::trial_space_type,
        typename FormType::trial_space_type>;

/** @brief Space-compatible canonical operator-pair predicate. */
template<typename ExprT, typename FormType,
         localform_basis_operator TestOperator,
         localform_basis_operator TrialOperator>
inline constexpr bool can_lower_canonical_bilinear_localform_for_form_v =
    can_lower_canonical_bilinear_localform_v<ExprT, TestOperator,
                                             TrialOperator> &&
    FormType::test_space_type::nSpaces == 1 &&
    FormType::trial_space_type::nSpaces == 1 &&
    std::is_same_v<
        typename canonical_bilinear_factors<ExprT>::test_space_type,
        typename FormType::test_space_type> &&
    std::is_same_v<
        typename canonical_bilinear_factors<ExprT>::trial_space_type,
        typename FormType::trial_space_type>;

/** @brief True for a vector value-inner-product mass form. */
template<typename ExprT, typename FormType>
inline constexpr bool can_lower_vector_mass_localform_for_form_v =
    can_lower_canonical_bilinear_localform_for_form_v<
        ExprT, FormType, localform_basis_operator::value,
        localform_basis_operator::value> &&
    ( !FormType::test_space_type::is_scalar ||
      !FormType::trial_space_type::is_scalar );

/** @brief True for an H(curl) curl-curl form. */
template<typename ExprT, typename FormType>
inline constexpr bool can_lower_curl_curl_localform_for_form_v =
    can_lower_canonical_bilinear_localform_for_form_v<
        ExprT, FormType, localform_basis_operator::curl,
        localform_basis_operator::curl> &&
    is_hcurl_conforming_v<
        typename FormType::test_space_type::reference_element_type> &&
    is_hcurl_conforming_v<
        typename FormType::trial_space_type::reference_element_type>;

/** @brief True for an H(div) divergence-divergence form. */
template<typename ExprT, typename FormType>
inline constexpr bool can_lower_div_div_localform_for_form_v =
    can_lower_canonical_bilinear_localform_for_form_v<
        ExprT, FormType, localform_basis_operator::divergence,
        localform_basis_operator::divergence> &&
    is_hdiv_conforming_v<
        typename FormType::test_space_type::reference_element_type> &&
    is_hdiv_conforming_v<
        typename FormType::trial_space_type::reference_element_type>;

/** @brief True for a scalar-test/H(div)-trial mixed divergence form. */
template<typename ExprT, typename FormType>
inline constexpr bool can_lower_value_div_localform_for_form_v =
    can_lower_canonical_bilinear_localform_for_form_v<
        ExprT, FormType, localform_basis_operator::value,
        localform_basis_operator::divergence> &&
    FormType::test_space_type::is_scalar &&
    is_hdiv_conforming_v<
        typename FormType::trial_space_type::reference_element_type>;

/** @brief True for an H(div)-test/scalar-trial mixed divergence form. */
template<typename ExprT, typename FormType>
inline constexpr bool can_lower_div_value_localform_for_form_v =
    can_lower_canonical_bilinear_localform_for_form_v<
        ExprT, FormType, localform_basis_operator::divergence,
        localform_basis_operator::value> &&
    is_hdiv_conforming_v<
        typename FormType::test_space_type::reference_element_type> &&
    FormType::trial_space_type::is_scalar;

} // namespace detail
} // namespace vf
} // namespace Feel

#undef FEELPP_LOCALFORM_GEOMETRY_SCALAR_LEAFS

#endif
