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
   \file localform_coeff_ir.hpp
   \brief Normalized scalar coefficient IR for lowered local forms.
 */
#ifndef FEELPP_VF_DETAIL_LOCALFORM_COEFF_IR_HPP
#define FEELPP_VF_DETAIL_LOCALFORM_COEFF_IR_HPP 1

#include <algorithm>

#include <feel/feelvf/detail/localform_base.hpp>

namespace Feel
{
namespace vf
{
namespace detail
{

struct localform_coeff_ir_marker {};

struct localform_coeff_constant_leaf_tag {};
struct localform_coeff_scalar_value_leaf_tag {};
struct localform_coeff_geometry_leaf_tag {};

struct localform_coeff_chi_op_tag {};
struct localform_coeff_add_op_tag {};
struct localform_coeff_mul_op_tag {};
struct localform_coeff_greater_op_tag {};

template<typename T>
struct is_localform_coeff_ir : std::is_base_of<localform_coeff_ir_marker, remove_cvref_t<T>> {};

template<typename T>
inline constexpr bool is_localform_coeff_ir_v = is_localform_coeff_ir<T>::value;

template<typename T>
using localform_coeff_evaluate_type_t = Eigen::Matrix<T, 1, 1>;

template<typename T>
inline auto
localform_scalar_from_evaluate( T const& eval )
{
    if constexpr ( std::is_arithmetic_v<T> )
        return eval;
    else
        return eval( 0, 0 );
}

template<typename ExprT, typename LeafTag>
class LocalformCoeffLeaf : public localform_coeff_ir_marker
{
public:
    using expression_type = remove_cvref_t<ExprT>;
    using leaf_tag = LeafTag;
    using value_type = typename expression_type::value_type;
    using evaluate_type = localform_coeff_evaluate_type_t<value_type>;
    static const size_type context = expression_type::context;
    static inline const bool is_terminal = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    static inline const bool has_test_basis = false;
    template<typename Func>
    static inline const bool has_trial_basis = false;

    explicit LocalformCoeffLeaf( expression_type expr )
        :
        M_expr( std::move( expr ) )
    {}

    size_type dynamicContext() const { return Feel::vf::dynamicContext( M_expr ); }
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }
    bool isPolynomial() const { return M_expr.isPolynomial(); }
    expression_type const& expression() const { return M_expr; }

    evaluate_type evaluate( bool p ) const
    {
        return evaluate_type::Constant( localform_scalar_from_evaluate( M_expr.evaluate( p ) ) );
    }

    template<typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_expr.applySymbolsExpr( se );
        using new_expr_type = remove_cvref_t<decltype( newExpr )>;
        return LocalformCoeffLeaf<new_expr_type, leaf_tag>( newExpr );
    }

    template<typename Geo_t, typename Basis_i_t = mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expr_type = expression_type;
        using inner_tensor_type = typename expr_type::template tensor<Geo_t>;
        using value_type = typename LocalformCoeffLeaf::value_type;
        using gmc_type = gmc_t<Geo_t>;
        using shape = Shape<gmc_type::nDim, Scalar, false, false>;
        using eval_matrix_type = Eigen::Matrix<value_type, shape::M, shape::N>;
        using eval_map_type = Eigen::Map<const eval_matrix_type>;

        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = inner_tensor_type::is_zero::value;
        };

        tensor( LocalformCoeffLeaf const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_expr( expr.expression(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}
        tensor( LocalformCoeffLeaf const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            M_expr( expr.expression(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}
        tensor( LocalformCoeffLeaf const& expr,
                Geo_t const& geom )
            :
            M_expr( expr.expression(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_expr.update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_expr.update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_expr.update( geom );
        }
        template<typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_expr.updateContext( ctx... );
        }

        value_type evalij( uint16_type i, uint16_type j ) const
        {
            return evalq( 0, 0, 0 );
        }
        eval_map_type evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2,
                            uint16_type q ) const
        {
            return M_expr.evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2,
                            uint16_type q,
                            mpl::int_<PatternContext> patternContext ) const
        {
            return M_expr.evalijq( i, j, c1, c2, q, patternContext );
        }
        value_type evaliq( uint16_type i,
                           uint16_type c1, uint16_type c2,
                           uint16_type q ) const
        {
            return M_expr.evalq( c1, c2, q );
        }
        eval_map_type evaliq( uint16_type i, uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }
        value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_expr.evalq( c1, c2, q );
        }
        eval_map_type evalq( uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }

        inner_tensor_type M_expr;
        mutable eval_matrix_type M_evalBuffer;
    };

private:
    expression_type M_expr;
};

template<typename ExprT>
using LocalformCoeffConstant = LocalformCoeffLeaf<ExprT, localform_coeff_constant_leaf_tag>;

template<typename ExprT>
using LocalformCoeffScalarValue = LocalformCoeffLeaf<ExprT, localform_coeff_scalar_value_leaf_tag>;

template<typename ExprT>
using LocalformCoeffGeometry = LocalformCoeffLeaf<ExprT, localform_coeff_geometry_leaf_tag>;

template<typename InnerExpr>
class LocalformCoeffChi : public localform_coeff_ir_marker
{
public:
    using inner_expr_type = remove_cvref_t<InnerExpr>;
    using value_type = typename inner_expr_type::value_type;
    using evaluate_type = localform_coeff_evaluate_type_t<value_type>;
    static const size_type context = inner_expr_type::context;
    static inline const bool is_terminal = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    static inline const bool has_test_basis = false;
    template<typename Func>
    static inline const bool has_trial_basis = false;

    explicit LocalformCoeffChi( inner_expr_type innerExpr )
        :
        M_innerExpr( std::move( innerExpr ) )
    {}

    size_type dynamicContext() const { return M_innerExpr.dynamicContext(); }
    uint16_type polynomialOrder() const { return M_innerExpr.polynomialOrder(); }
    bool isPolynomial() const { return M_innerExpr.isPolynomial() && M_innerExpr.polynomialOrder() == 0; }
    inner_expr_type const& expression() const { return M_innerExpr; }

    evaluate_type evaluate( bool p ) const
    {
        auto value = localform_scalar_from_evaluate( M_innerExpr.evaluate( p ) );
        return evaluate_type::Constant( value != value_type( 0 ) ? value_type( 1 ) : value_type( 0 ) );
    }

    template<typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_innerExpr.applySymbolsExpr( se );
        using new_expr_type = remove_cvref_t<decltype( newExpr )>;
        return LocalformCoeffChi<new_expr_type>( newExpr );
    }

    template<typename Geo_t, typename Basis_i_t = mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using inner_tensor_type = typename inner_expr_type::template tensor<Geo_t>;
        using value_type = typename LocalformCoeffChi::value_type;
        using gmc_type = gmc_t<Geo_t>;
        using shape = Shape<gmc_type::nDim, Scalar, false, false>;
        using eval_matrix_type = Eigen::Matrix<value_type, shape::M, shape::N>;
        using eval_map_type = Eigen::Map<const eval_matrix_type>;

        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = inner_tensor_type::is_zero::value;
        };

        tensor( LocalformCoeffChi const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_innerExpr( expr.expression(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}
        tensor( LocalformCoeffChi const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            M_innerExpr( expr.expression(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}
        tensor( LocalformCoeffChi const& expr,
                Geo_t const& geom )
            :
            M_innerExpr( expr.expression(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_innerExpr.update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_innerExpr.update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_innerExpr.update( geom );
        }
        template<typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_innerExpr.updateContext( ctx... );
        }

        value_type evalij( uint16_type i, uint16_type j ) const
        {
            return evalq( 0, 0, 0 );
        }
        eval_map_type evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2,
                            uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2,
                            uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1, c2, q );
        }
        value_type evaliq( uint16_type i,
                           uint16_type c1, uint16_type c2,
                           uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        eval_map_type evaliq( uint16_type i, uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }
        value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_innerExpr.evalq( c1, c2, q ) != value_type( 0 ) ? value_type( 1 ) : value_type( 0 );
        }
        eval_map_type evalq( uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }

        inner_tensor_type M_innerExpr;
        mutable eval_matrix_type M_evalBuffer;
    };

private:
    inner_expr_type M_innerExpr;
};

template<typename BinaryOpTag, typename LeftValueType, typename RightValueType>
struct LocalformCoeffBinaryOp;

template<typename LeftValueType, typename RightValueType>
struct LocalformCoeffBinaryOp<localform_coeff_add_op_tag, LeftValueType, RightValueType>
{
    using value_type = std::common_type_t<LeftValueType, RightValueType>;

    static value_type apply( LeftValueType const& left, RightValueType const& right )
    {
        return left + right;
    }
    static uint16_type polynomialOrder( uint16_type left, uint16_type right )
    {
        return std::max( left, right );
    }
    static bool isPolynomial( bool leftIsPolynomial, bool rightIsPolynomial )
    {
        return leftIsPolynomial && rightIsPolynomial;
    }
};

template<typename LeftValueType, typename RightValueType>
struct LocalformCoeffBinaryOp<localform_coeff_mul_op_tag, LeftValueType, RightValueType>
{
    using value_type = std::common_type_t<LeftValueType, RightValueType>;

    static value_type apply( LeftValueType const& left, RightValueType const& right )
    {
        return left * right;
    }
    static uint16_type polynomialOrder( uint16_type left, uint16_type right )
    {
        return left + right;
    }
    static bool isPolynomial( bool leftIsPolynomial, bool rightIsPolynomial )
    {
        return leftIsPolynomial && rightIsPolynomial;
    }
};

template<typename LeftValueType, typename RightValueType>
struct LocalformCoeffBinaryOp<localform_coeff_greater_op_tag, LeftValueType, RightValueType>
{
    using value_type = std::common_type_t<LeftValueType, RightValueType>;

    static value_type apply( LeftValueType const& left, RightValueType const& right )
    {
        return left > right ? value_type( 1 ) : value_type( 0 );
    }
    static uint16_type polynomialOrder( uint16_type /*left*/, uint16_type /*right*/ )
    {
        return 0;
    }
    static bool isPolynomial( bool /*leftIsPolynomial*/, bool /*rightIsPolynomial*/ )
    {
        return false;
    }
};

template<typename LeftExpr, typename RightExpr, typename BinaryOpTag>
class LocalformCoeffBinary : public localform_coeff_ir_marker
{
public:
    using left_expr_type = remove_cvref_t<LeftExpr>;
    using right_expr_type = remove_cvref_t<RightExpr>;
    using op_type = LocalformCoeffBinaryOp<BinaryOpTag, typename left_expr_type::value_type, typename right_expr_type::value_type>;
    using value_type = typename op_type::value_type;
    using evaluate_type = localform_coeff_evaluate_type_t<value_type>;
    static const size_type context = left_expr_type::context | right_expr_type::context;
    static inline const bool is_terminal = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    template<typename Func>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };
    template<typename Func>
    static inline const bool has_test_basis = false;
    template<typename Func>
    static inline const bool has_trial_basis = false;

    LocalformCoeffBinary( left_expr_type leftExpr, right_expr_type rightExpr )
        :
        M_leftExpr( std::move( leftExpr ) ),
        M_rightExpr( std::move( rightExpr ) )
    {}

    size_type dynamicContext() const { return M_leftExpr.dynamicContext() | M_rightExpr.dynamicContext(); }
    uint16_type polynomialOrder() const
    {
        return op_type::polynomialOrder( M_leftExpr.polynomialOrder(), M_rightExpr.polynomialOrder() );
    }
    bool isPolynomial() const
    {
        return op_type::isPolynomial( M_leftExpr.isPolynomial(), M_rightExpr.isPolynomial() );
    }
    left_expr_type const& left() const { return M_leftExpr; }
    right_expr_type const& right() const { return M_rightExpr; }

    evaluate_type evaluate( bool p ) const
    {
        auto left = localform_scalar_from_evaluate( M_leftExpr.evaluate( p ) );
        auto right = localform_scalar_from_evaluate( M_rightExpr.evaluate( p ) );
        return evaluate_type::Constant( op_type::apply( left, right ) );
    }

    template<typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newLeftExpr = M_leftExpr.applySymbolsExpr( se );
        auto newRightExpr = M_rightExpr.applySymbolsExpr( se );
        using new_left_expr_type = remove_cvref_t<decltype( newLeftExpr )>;
        using new_right_expr_type = remove_cvref_t<decltype( newRightExpr )>;
        return LocalformCoeffBinary<new_left_expr_type, new_right_expr_type, BinaryOpTag>( newLeftExpr, newRightExpr );
    }

    template<typename Geo_t, typename Basis_i_t = mpl::void_, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using left_tensor_type = typename left_expr_type::template tensor<Geo_t>;
        using right_tensor_type = typename right_expr_type::template tensor<Geo_t>;
        using value_type = typename LocalformCoeffBinary::value_type;
        using gmc_type = gmc_t<Geo_t>;
        using shape = Shape<gmc_type::nDim, Scalar, false, false>;
        using eval_matrix_type = Eigen::Matrix<value_type, shape::M, shape::N>;
        using eval_map_type = Eigen::Map<const eval_matrix_type>;

        template<typename Indq, typename Indi, typename Indj>
        struct expr
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( LocalformCoeffBinary const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/,
                Basis_j_t const& /*feu*/ )
            :
            M_leftExpr( expr.left(), geom ),
            M_rightExpr( expr.right(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}
        tensor( LocalformCoeffBinary const& expr,
                Geo_t const& geom,
                Basis_i_t const& /*fev*/ )
            :
            M_leftExpr( expr.left(), geom ),
            M_rightExpr( expr.right(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}
        tensor( LocalformCoeffBinary const& expr,
                Geo_t const& geom )
            :
            M_leftExpr( expr.left(), geom ),
            M_rightExpr( expr.right(), geom ),
            M_evalBuffer( eval_matrix_type::Zero() )
        {}

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_leftExpr.update( geom );
            M_rightExpr.update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_leftExpr.update( geom );
            M_rightExpr.update( geom );
        }
        void update( Geo_t const& geom )
        {
            M_leftExpr.update( geom );
            M_rightExpr.update( geom );
        }
        template<typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_leftExpr.updateContext( ctx... );
            M_rightExpr.updateContext( ctx... );
        }

        value_type evalij( uint16_type i, uint16_type j ) const
        {
            return evalq( 0, 0, 0 );
        }
        eval_map_type evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2,
                            uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        template<int PatternContext>
        value_type evalijq( uint16_type i, uint16_type j,
                            uint16_type c1, uint16_type c2,
                            uint16_type q,
                            mpl::int_<PatternContext> ) const
        {
            return evalijq( i, j, c1, c2, q );
        }
        value_type evaliq( uint16_type i,
                           uint16_type c1, uint16_type c2,
                           uint16_type q ) const
        {
            return evalq( c1, c2, q );
        }
        eval_map_type evaliq( uint16_type i, uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }
        value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return op_type::apply( M_leftExpr.evalq( c1, c2, q ), M_rightExpr.evalq( c1, c2, q ) );
        }
        eval_map_type evalq( uint16_type q ) const
        {
            M_evalBuffer( 0, 0 ) = evalq( 0, 0, q );
            return eval_map_type( M_evalBuffer.data() );
        }

        left_tensor_type M_leftExpr;
        right_tensor_type M_rightExpr;
        mutable eval_matrix_type M_evalBuffer;
    };

private:
    left_expr_type M_leftExpr;
    right_expr_type M_rightExpr;
};

template<typename LeftExpr, typename RightExpr>
using LocalformCoeffAdd = LocalformCoeffBinary<LeftExpr, RightExpr, localform_coeff_add_op_tag>;

template<typename LeftExpr, typename RightExpr>
using LocalformCoeffMul = LocalformCoeffBinary<LeftExpr, RightExpr, localform_coeff_mul_op_tag>;

template<typename LeftExpr, typename RightExpr>
using LocalformCoeffGreater = LocalformCoeffBinary<LeftExpr, RightExpr, localform_coeff_greater_op_tag>;

} // namespace detail
} // namespace vf
} // namespace Feel

#endif
