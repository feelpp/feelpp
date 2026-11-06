/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme
       Date: 2026-03-25

  Copyright (C) 2026 Feel++ Consortium

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file tensorbasis.hpp
   \author Christophe Prud'homme
   \date 2026-03-25
 */
#ifndef FEELPP_VF_TENSORBASIS_HPP
#define FEELPP_VF_TENSORBASIS_HPP 1

#include <concepts>
#include <numbers>
#include <type_traits>

#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/matvec.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/ones.hpp>
#include <feel/feelvf/symm.hpp>
#include <feel/feelvf/trans.hpp>

namespace Feel
{
namespace vf
{

namespace detail
{

enum class TensorBasisKind
{
    Canonical,
    Symmetric,
    Mandel
};

template <int Dim, int I, int K, TensorBasisKind Kind>
class TensorBasis
{
public:
    static_assert( Dim > 0, "TensorBasis expects a strictly positive dimension" );
    static_assert( I >= 0 && I < Dim, "TensorBasis row index is out of bounds" );
    static_assert( K >= 0 && K < Dim, "TensorBasis column index is out of bounds" );

    static constexpr size_type context = 0;
    static inline const bool is_terminal = true;

    template <typename Func>
    struct HasTestFunction
    {
        static inline const bool result = false;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = false;
    };

    template <typename Func>
    static inline const bool has_test_basis = false;
    template <typename Func>
    static inline const bool has_trial_basis = false;
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    using this_type = TensorBasis<Dim, I, K, Kind>;
    using value_type = double;
    using evaluate_type = Eigen::Matrix<value_type, Dim, Dim>;

    static constexpr value_type mandel_scale = value_type( 1 )/std::numbers::sqrt2_v<value_type>;

    TensorBasis() = default;
    TensorBasis( TensorBasis const& ) = default;
    TensorBasis& operator=( TensorBasis const& ) = default;

    template <typename... TheExpr>
    struct Lambda
    {
        using type = this_type;
    };

    template <typename... TheExpr>
    [[nodiscard]] this_type operator()( TheExpr... ) const
    {
        return {};
    }

    constexpr uint16_type polynomialOrder() const { return 0; }
    constexpr bool isPolynomial() const { return true; }

    evaluate_type evaluate( bool ) const
    {
        evaluate_type result = evaluate_type::Zero();
        for ( int c1 = 0; c1 < Dim; ++c1 )
            for ( int c2 = 0; c2 < Dim; ++c2 )
                result( c1, c2 ) = basisValue( c1, c2 );
        return result;
    }

    template <typename SymbolsExprType>
    [[nodiscard]] this_type applySymbolsExpr( SymbolsExprType const& ) const
    {
        return {};
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const&, WorldComm const&, std::string const&, TheSymbolExprType const& ) const
    {
        return zero<Dim, Dim>();
    }

    [[nodiscard]] static constexpr value_type basisValue( uint16_type c1, uint16_type c2 )
    {
        if constexpr ( Kind == TensorBasisKind::Canonical )
        {
            return ( c1 == I && c2 == K ) ? value_type( 1 ) : value_type( 0 );
        }
        else if constexpr ( I == K )
        {
            return ( c1 == I && c2 == K ) ? value_type( 1 ) : value_type( 0 );
        }
        else if constexpr ( Kind == TensorBasisKind::Symmetric )
        {
            return ( ( c1 == I && c2 == K ) || ( c1 == K && c2 == I ) ) ? value_type( 1 ) : value_type( 0 );
        }
        else
        {
            return ( ( c1 == I && c2 == K ) || ( c1 == K && c2 == I ) ) ? mandel_scale : value_type( 0 );
        }
    }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using expression_type = this_type;
        using key_type = key_t<Geo_t>;
        using gmc_type = gmc_t<Geo_t>;
        using shape = ShapeGeneric<gmc_type::nDim, Dim, Dim>;
        using value_type = typename expression_type::value_type;

        template <class Args>
        struct sig
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = false;
        };

        tensor( expression_type const&, Geo_t const&, Basis_i_t const&, Basis_j_t const& )
        {
        }

        tensor( expression_type const&, Geo_t const&, Basis_i_t const& )
        {
        }

        tensor( expression_type const&, Geo_t const& )
        {
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&, expression_type const&,
                Geo_t const&, TheArgsType const&... )
        {
        }

        void update( Geo_t const&, Basis_i_t const&, Basis_j_t const& )
        {
        }

        void update( Geo_t const&, Basis_i_t const& )
        {
        }

        void update( Geo_t const& )
        {
        }

        template <typename... CTX>
        void updateContext( CTX const&... )
        {
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type, TheExprExpandedType const&, TupleTensorSymbolsExprType&, Geo_t const&,
                     TheArgsType const&... )
        {
        }

        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type c2, uint16_type ) const
        {
            return expression_type::basisValue( c1, c2 );
        }

        template <int PatternContext>
        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type, uint16_type, uint16_type c1, uint16_type c2, uint16_type, mpl::int_<PatternContext> ) const
        {
            return expression_type::basisValue( c1, c2 );
        }

        FEELPP_STRONG_INLINE value_type
        evaliq( uint16_type, uint16_type c1, uint16_type c2, uint16_type ) const
        {
            return expression_type::basisValue( c1, c2 );
        }

        FEELPP_STRONG_INLINE value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type ) const
        {
            return expression_type::basisValue( c1, c2 );
        }
    };
};

template <int Dim, int I, int K, TensorBasisKind Kind>
using tensor_basis_expr_t = Expr<TensorBasis<Dim, I, K, Kind>>;

template <int Dim, int I, int K, TensorBasisKind Kind>
[[nodiscard]] inline tensor_basis_expr_t<Dim, I, K, Kind>
makeTensorBasisExpr()
{
    using basis_type = TensorBasis<Dim, I, K, Kind>;
    return tensor_basis_expr_t<Dim, I, K, Kind>( basis_type{} );
}

template <typename T>
struct is_tensor_basis_expression : std::false_type
{
};

template <int Dim, int I, int K, TensorBasisKind Kind>
struct is_tensor_basis_expression<tensor_basis_expr_t<Dim, I, K, Kind>> : std::true_type
{
};

template <typename T>
concept TensorBasisExpression = is_tensor_basis_expression<std::remove_cvref_t<T>>::value;

template <typename ExprT, int Dim, int I, int K>
class LeftMultiplyByDelta : public ExprDynamicBase
{
public:
    using super = ExprDynamicBase;
    static const size_type context = ExprT::context;
    static inline const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static inline const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    template <typename Func>
    static inline const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template <typename Func>
    static inline const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    using expression_type = ExprT;
    using value_type = typename expression_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type,
                                        Dim,
                                        expression_type::evaluate_type::ColsAtCompileTime>;
    using this_type = LeftMultiplyByDelta<ExprT, Dim, I, K>;

    explicit LeftMultiplyByDelta( expression_type const& expr )
        :
        super( Feel::vf::dynamicContext( expr ) ),
        M_expr( expr )
    {
    }

    template <typename... TheExpr>
    struct Lambda
    {
        using type = LeftMultiplyByDelta<typename expression_type::template Lambda<TheExpr...>::type, Dim, I, K>;
    };

    template <typename... TheExpr>
    [[nodiscard]] typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) const
    {
        return typename Lambda<TheExpr...>::type( M_expr( e... ) );
    }

    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    evaluate_type evaluate( bool p ) const
    {
        auto source = M_expr.evaluate( p ).template cast<value_type>();
        evaluate_type result;
        if constexpr ( evaluate_type::RowsAtCompileTime == Eigen::Dynamic ||
                       evaluate_type::ColsAtCompileTime == Eigen::Dynamic )
            result = evaluate_type::Zero( Dim, source.cols() );
        else
            result = evaluate_type::Zero();
        result.row( I ) = source.row( K );
        return result;
    }

    void setParameterValues( std::map<std::string, value_type> const& mp )
    {
        M_expr.setParameterValues( mp );
    }

    void updateParameterValues( std::map<std::string, double>& pv ) const
    {
        M_expr.updateParameterValues( pv );
    }

    template <typename SymbolsExprType>
    [[nodiscard]] auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_expr.applySymbolsExpr( se );
        using new_expr_type = std::decay_t<decltype( newExpr )>;
        return LeftMultiplyByDelta<new_expr_type, Dim, I, K>( newExpr );
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return M_expr.hasSymbolDependency( symb, se );
    }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string, std::set<std::string>>& res,
                           TheSymbolExprType const& se ) const
    {
        M_expr.dependentSymbols( symb, res, se );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto diffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
        using diff_expr_type = std::decay_t<decltype( diffExpr )>;
        return LeftMultiplyByDelta<diff_expr_type, Dim, I, K>( diffExpr );
    }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using tensor_expr_type = typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using value_type = typename tensor_expr_type::value_type;
        using expr_shape = typename tensor_expr_type::shape;
        using shape = ShapeGeneric<expr_shape::nDim, Dim, expr_shape::N>;

        BOOST_MPL_ASSERT_MSG( ( expr_shape::M == Dim ),
                              LEFT_MULTIPLICATION_BY_DELTA_REQUIRES_A_COMPATIBLE_ROW_DIMENSION,
                              ( mpl::int_<expr_shape::M>, mpl::int_<Dim> ) );

        template <class Args>
        struct sig
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                this_type const& expr, Geo_t const& geom, TheArgsType const&... theInitArgs )
            :
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
        {
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }

        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }

        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }

        template <typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_tensor_expr.updateContext( ctx... );
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                     Geo_t const& geom, TheArgsType const&... theUpdateArgs )
        {
            M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
        }

        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return ( c1 == I ) ? M_tensor_expr.evalijq( i, j, K, c2, q ) : value_type( 0 );
        }

        template <int PatternContext>
        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return ( c1 == I ) ? M_tensor_expr.evalijq( i, j, K, c2, q, mpl::int_<PatternContext>() ) : value_type( 0 );
        }

        FEELPP_STRONG_INLINE value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return ( c1 == I ) ? M_tensor_expr.evaliq( i, K, c2, q ) : value_type( 0 );
        }

        FEELPP_STRONG_INLINE value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return ( c1 == I ) ? M_tensor_expr.evalq( K, c2, q ) : value_type( 0 );
        }

    private:
        tensor_expr_type M_tensor_expr;
    };

private:
    mutable expression_type M_expr;
};

template <typename ExprT, int Dim, int I, int K>
class RightMultiplyByDelta : public ExprDynamicBase
{
public:
    using super = ExprDynamicBase;
    static const size_type context = ExprT::context;
    static inline const bool is_terminal = false;

    template <typename Func>
    struct HasTestFunction
    {
        static inline const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static inline const bool result = ExprT::template HasTrialFunction<Func>::result;
    };

    template <typename Func>
    static inline const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template <typename Func>
    static inline const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    using expression_type = ExprT;
    using value_type = typename expression_type::value_type;
    using evaluate_type = Eigen::Matrix<value_type,
                                        expression_type::evaluate_type::RowsAtCompileTime,
                                        Dim>;
    using this_type = RightMultiplyByDelta<ExprT, Dim, I, K>;

    explicit RightMultiplyByDelta( expression_type const& expr )
        :
        super( Feel::vf::dynamicContext( expr ) ),
        M_expr( expr )
    {
    }

    template <typename... TheExpr>
    struct Lambda
    {
        using type = RightMultiplyByDelta<typename expression_type::template Lambda<TheExpr...>::type, Dim, I, K>;
    };

    template <typename... TheExpr>
    [[nodiscard]] typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) const
    {
        return typename Lambda<TheExpr...>::type( M_expr( e... ) );
    }

    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    expression_type const& expression() const
    {
        return M_expr;
    }

    evaluate_type evaluate( bool p ) const
    {
        auto source = M_expr.evaluate( p ).template cast<value_type>();
        evaluate_type result;
        if constexpr ( evaluate_type::RowsAtCompileTime == Eigen::Dynamic ||
                       evaluate_type::ColsAtCompileTime == Eigen::Dynamic )
            result = evaluate_type::Zero( source.rows(), Dim );
        else
            result = evaluate_type::Zero();
        result.col( K ) = source.col( I );
        return result;
    }

    void setParameterValues( std::map<std::string, value_type> const& mp )
    {
        M_expr.setParameterValues( mp );
    }

    void updateParameterValues( std::map<std::string, double>& pv ) const
    {
        M_expr.updateParameterValues( pv );
    }

    template <typename SymbolsExprType>
    [[nodiscard]] auto applySymbolsExpr( SymbolsExprType const& se ) const
    {
        auto newExpr = M_expr.applySymbolsExpr( se );
        using new_expr_type = std::decay_t<decltype( newExpr )>;
        return RightMultiplyByDelta<new_expr_type, Dim, I, K>( newExpr );
    }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
    {
        return M_expr.hasSymbolDependency( symb, se );
    }

    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string, std::set<std::string>>& res,
                           TheSymbolExprType const& se ) const
    {
        M_expr.dependentSymbols( symb, res, se );
    }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
    {
        auto diffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
        using diff_expr_type = std::decay_t<decltype( diffExpr )>;
        return RightMultiplyByDelta<diff_expr_type, Dim, I, K>( diffExpr );
    }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        using tensor_expr_type = typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>;
        using value_type = typename tensor_expr_type::value_type;
        using expr_shape = typename tensor_expr_type::shape;
        using shape = ShapeGeneric<expr_shape::nDim, expr_shape::M, Dim>;

        BOOST_MPL_ASSERT_MSG( ( expr_shape::N == Dim ),
                              RIGHT_MULTIPLICATION_BY_DELTA_REQUIRES_A_COMPATIBLE_COLUMN_DIMENSION,
                              ( mpl::int_<expr_shape::N>, mpl::int_<Dim> ) );

        template <class Args>
        struct sig
        {
            using type = value_type;
        };

        struct is_zero
        {
            static inline const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_tensor_expr( expr.expression(), geom, fev, feu )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            M_tensor_expr( expr.expression(), geom, fev )
        {
        }

        tensor( this_type const& expr, Geo_t const& geom )
            :
            M_tensor_expr( expr.expression(), geom )
        {
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                this_type const& expr, Geo_t const& geom, TheArgsType const&... theInitArgs )
            :
            M_tensor_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
        {
        }

        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensor_expr.update( geom, fev, feu );
        }

        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensor_expr.update( geom, fev );
        }

        void update( Geo_t const& geom )
        {
            M_tensor_expr.update( geom );
        }

        template <typename... CTX>
        void updateContext( CTX const&... ctx )
        {
            M_tensor_expr.updateContext( ctx... );
        }

        template <typename TheExprExpandedType, typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType& ttse,
                     Geo_t const& geom, TheArgsType const&... theUpdateArgs )
        {
            M_tensor_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
        }

        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return ( c2 == K ) ? M_tensor_expr.evalijq( i, j, c1, I, q ) : value_type( 0 );
        }

        template <int PatternContext>
        FEELPP_STRONG_INLINE value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return ( c2 == K ) ? M_tensor_expr.evalijq( i, j, c1, I, q, mpl::int_<PatternContext>() ) : value_type( 0 );
        }

        FEELPP_STRONG_INLINE value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return ( c2 == K ) ? M_tensor_expr.evaliq( i, c1, I, q ) : value_type( 0 );
        }

        FEELPP_STRONG_INLINE value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return ( c2 == K ) ? M_tensor_expr.evalq( c1, I, q ) : value_type( 0 );
        }

    private:
        tensor_expr_type M_tensor_expr;
    };

private:
    mutable expression_type M_expr;
};

template <typename ExprT, int Dim, int I, int K>
[[nodiscard]] inline auto
makeLeftMultiplyByDeltaExpr( ExprT const& expr )
{
    using expr_type = LeftMultiplyByDelta<std::decay_t<ExprT>, Dim, I, K>;
    return Expr<expr_type>( expr_type( expr ) );
}

template <typename ExprT, int Dim, int I, int K>
[[nodiscard]] inline auto
makeRightMultiplyByDeltaExpr( ExprT const& expr )
{
    using expr_type = RightMultiplyByDelta<std::decay_t<ExprT>, Dim, I, K>;
    return Expr<expr_type>( expr_type( expr ) );
}

template <int Dim, int I, int K, TensorBasisKind Kind, typename ExprT>
[[nodiscard]] inline auto
tensorBasisInnerExpr( ExprT const& expr )
{
    if constexpr ( Kind == TensorBasisKind::Canonical || I == K )
    {
        return component<I, K>( expr );
    }
    else if constexpr ( Kind == TensorBasisKind::Symmetric )
    {
        return component<I, K>( expr ) + component<K, I>( expr );
    }
    else
    {
        return cst( TensorBasis<Dim, I, K, Kind>::mandel_scale )*( component<I, K>( expr ) + component<K, I>( expr ) );
    }
}

template <int Dim1, int I1, int K1, TensorBasisKind Kind1, int Dim2, int I2, int K2, TensorBasisKind Kind2>
[[nodiscard]] consteval double
tensorBasisInnerValue()
{
    static_assert( Dim1 == Dim2, "inner() between tensor bases expects matching dimensions" );

    double result = 0.0;
    for ( int c1 = 0; c1 < Dim1; ++c1 )
        for ( int c2 = 0; c2 < Dim1; ++c2 )
            result += TensorBasis<Dim1, I1, K1, Kind1>::basisValue( c1, c2 ) *
                      TensorBasis<Dim2, I2, K2, Kind2>::basisValue( c1, c2 );
    return result;
}

template <int Dim, TensorBasisKind Kind>
[[nodiscard]] inline double
dynamicTensorBasisValue( int i, int k, int c1, int c2 )
{
    CHECK( i >= 0 && i < Dim ) << "runtime tensor basis row index out of bounds: " << i;
    CHECK( k >= 0 && k < Dim ) << "runtime tensor basis column index out of bounds: " << k;

    if constexpr ( Kind == TensorBasisKind::Canonical )
    {
        return ( c1 == i && c2 == k ) ? 1.0 : 0.0;
    }
    else if ( i == k )
    {
        return ( c1 == i && c2 == k ) ? 1.0 : 0.0;
    }
    else if constexpr ( Kind == TensorBasisKind::Symmetric )
    {
        return ( ( c1 == i && c2 == k ) || ( c1 == k && c2 == i ) ) ? 1.0 : 0.0;
    }
    else
    {
        return ( ( c1 == i && c2 == k ) || ( c1 == k && c2 == i ) ) ?
            1.0/std::numbers::sqrt2_v<double> : 0.0;
    }
}

template <int Dim, TensorBasisKind Kind, std::size_t... FlatIndices>
[[nodiscard]] inline auto
makeDynamicTensorBasisExpr( int i, int k, std::index_sequence<FlatIndices...> )
{
    return mat<Dim, Dim>( cst( dynamicTensorBasisValue<Dim, Kind>( i, k,
                                                                   static_cast<int>( FlatIndices/Dim ),
                                                                   static_cast<int>( FlatIndices%Dim ) ) )... );
}

} // namespace detail

template <int Dim, int I, int K>
[[nodiscard]] inline auto
delta()
{
    return detail::makeTensorBasisExpr<Dim, I, K, detail::TensorBasisKind::Canonical>();
}

template <int Dim>
[[nodiscard]] inline auto
delta( int i, int k )
{
    return detail::makeDynamicTensorBasisExpr<Dim, detail::TensorBasisKind::Canonical>( i, k,
                                                                                         std::make_index_sequence<Dim*Dim>{} );
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
symm_delta()
{
    return detail::makeTensorBasisExpr<Dim, I, K, detail::TensorBasisKind::Symmetric>();
}

template <int Dim>
[[nodiscard]] inline auto
symm_delta( int i, int k )
{
    return detail::makeDynamicTensorBasisExpr<Dim, detail::TensorBasisKind::Symmetric>( i, k,
                                                                                         std::make_index_sequence<Dim*Dim>{} );
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
mandel_delta()
{
    return detail::makeTensorBasisExpr<Dim, I, K, detail::TensorBasisKind::Mandel>();
}

template <int Dim>
[[nodiscard]] inline auto
mandel_delta( int i, int k )
{
    return detail::makeDynamicTensorBasisExpr<Dim, detail::TensorBasisKind::Mandel>( i, k,
                                                                                      std::make_index_sequence<Dim*Dim>{} );
}

template <int Dim, int I, int K, detail::TensorBasisKind Kind>
[[nodiscard]] inline auto
trans( detail::tensor_basis_expr_t<Dim, I, K, Kind> const& )
{
    return detail::makeTensorBasisExpr<Dim, K, I, Kind>();
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
sym( detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Canonical> const& )
{
    if constexpr ( I == K )
        return delta<Dim, I, K>();
    else
        return cst( 0.5 )*symm_delta<Dim, I, K>();
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
sym( detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Symmetric> const& basis )
{
    return basis;
}

template <int Dim, int I, int K>
[[nodiscard]] inline auto
sym( detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Mandel> const& basis )
{
    return basis;
}

template <int Dim, int I, int K, typename ExprT>
    requires ExpressionRowsCompatible<Expr<ExprT>, Dim>
[[nodiscard]] inline auto
operator*( detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Canonical> const&, Expr<ExprT> const& expr )
{
    return detail::makeLeftMultiplyByDeltaExpr<Expr<ExprT>, Dim, I, K>( expr );
}

template <typename ExprT, int Dim, int I, int K>
    requires ExpressionColsCompatible<Expr<ExprT>, Dim>
[[nodiscard]] inline auto
operator*( Expr<ExprT> const& expr, detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Canonical> const& )
{
    return detail::makeRightMultiplyByDeltaExpr<Expr<ExprT>, Dim, I, K>( expr );
}

template <int Dim, int I, int K, typename ExprT>
    requires ExpressionRowsCompatible<Expr<ExprT>, Dim>
[[nodiscard]] inline auto
operator*( detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Symmetric> const&, Expr<ExprT> const& expr )
{
    if constexpr ( I == K )
        return delta<Dim, I, K>()*expr;
    else
        return delta<Dim, I, K>()*expr + delta<Dim, K, I>()*expr;
}

template <typename ExprT, int Dim, int I, int K>
    requires ExpressionColsCompatible<Expr<ExprT>, Dim>
[[nodiscard]] inline auto
operator*( Expr<ExprT> const& expr, detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Symmetric> const& )
{
    if constexpr ( I == K )
        return expr*delta<Dim, I, K>();
    else
        return expr*delta<Dim, I, K>() + expr*delta<Dim, K, I>();
}

template <int Dim, int I, int K, typename ExprT>
    requires ExpressionRowsCompatible<Expr<ExprT>, Dim>
[[nodiscard]] inline auto
operator*( detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Mandel> const&, Expr<ExprT> const& expr )
{
    if constexpr ( I == K )
        return delta<Dim, I, K>()*expr;
    else
        return cst( detail::TensorBasis<Dim, I, K, detail::TensorBasisKind::Mandel>::mandel_scale )*
               ( delta<Dim, I, K>()*expr + delta<Dim, K, I>()*expr );
}

template <typename ExprT, int Dim, int I, int K>
    requires ExpressionColsCompatible<Expr<ExprT>, Dim>
[[nodiscard]] inline auto
operator*( Expr<ExprT> const& expr, detail::tensor_basis_expr_t<Dim, I, K, detail::TensorBasisKind::Mandel> const& )
{
    if constexpr ( I == K )
        return expr*delta<Dim, I, K>();
    else
        return cst( detail::TensorBasis<Dim, I, K, detail::TensorBasisKind::Mandel>::mandel_scale )*
               ( expr*delta<Dim, I, K>() + expr*delta<Dim, K, I>() );
}

template <typename ExprT, int Dim, int I, int K, detail::TensorBasisKind Kind>
    requires MatrixComponentAccessible<ExprT>
[[nodiscard]] inline auto
inner( ExprT const& expr, detail::tensor_basis_expr_t<Dim, I, K, Kind> const& )
{
    return detail::tensorBasisInnerExpr<Dim, I, K, Kind>( expr );
}

template <int Dim, int I, int K, detail::TensorBasisKind Kind, typename ExprT>
    requires MatrixComponentAccessible<ExprT>
[[nodiscard]] inline auto
inner( detail::tensor_basis_expr_t<Dim, I, K, Kind> const&, ExprT const& expr )
{
    return detail::tensorBasisInnerExpr<Dim, I, K, Kind>( expr );
}

template <int Dim1, int I1, int K1, detail::TensorBasisKind Kind1,
          int Dim2, int I2, int K2, detail::TensorBasisKind Kind2>
[[nodiscard]] inline auto
inner( detail::tensor_basis_expr_t<Dim1, I1, K1, Kind1> const&,
       detail::tensor_basis_expr_t<Dim2, I2, K2, Kind2> const& )
{
    return cst( detail::tensorBasisInnerValue<Dim1, I1, K1, Kind1, Dim2, I2, K2, Kind2>() );
}

} // namespace vf
} // namespace Feel

#endif /* FEELPP_VF_TENSORBASIS_HPP */
