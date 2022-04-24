/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-05-13

  Copyright (C) 2014-2016 Feel++ Consortium

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
   \file unary.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-05-13
 */
#ifndef FEELPP_VF_UNARY_HPP
#define FEELPP_VF_UNARY_HPP 1

namespace Feel
{
namespace vf
{
/**
   \class UnaryPlus
   \brief handler for unary plus expression
*/
template < class ExprT >
class UnaryPlus : public ExprDynamicBase
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    typedef ExprT expression_type;
    typedef typename ExprT::value_type value_type;
    using evaluate_type = typename expression_type::evaluate_type;
    typedef UnaryPlus<ExprT> this_type;

    UnaryPlus( const ExprT& expr )
        :
        M_expr( expr )
    {
    }
    UnaryPlus( UnaryPlus const& e )
        :
        M_expr( e.M_expr )
    {
    }

    template<typename... TheExpr>
    struct Lambda
    {
        typedef UnaryPlus<typename ExprT::template Lambda<TheExpr...>::type> type;
    };

    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_expr( e... ) ); }

    size_type dynamicContext() const { return vf::dynamicContext( M_expr ); } 

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            return M_expr.evaluate(p,worldcomm);
        }

    expression_type const& expression() const
    {
        return M_expr;
    }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            M_expr.setParameterValues( mp );
        }
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            M_expr.updateParameterValues( pv );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newExpr = M_expr.applySymbolsExpr( se );
            using new_expr_type = std::decay_t<decltype(newExpr)>;
            return UnaryPlus<new_expr_type>( newExpr );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr.hasSymbolDependency( symb, se );
        }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_expr.dependentSymbols( symb, res, se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            auto theDiffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
            using new_expr_type = std::decay_t<decltype(theDiffExpr)>;
            return UnaryPlus<new_expr_type>( theDiffExpr );
        }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape shape;

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_t_expr( expr.expression(), geom, fev, feu )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_t_expr( expr.expression(), geom, fev )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_t_expr( expr.expression(), geom )
        {}
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_t_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
            {}

        void update( Geo_t const& geom, Basis_i_t const& fev , Basis_j_t const& feu )
        {
            M_t_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev  )
        {
            M_t_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_t_expr.update( geom );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_t_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
            }
        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            M_t_expr.updateContext( ctx... );
        }

        value_type
        evalij( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2 ) const
        {
            return M_t_expr.evalij( i, j, c1, c2 );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_t_expr.evalijq( i, j, c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return M_t_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_t_expr.evaliq( i, c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_t_expr.evalq( c1, c2, q );
        }
        tensor_expr_type M_t_expr;
    };
protected:
    UnaryPlus() {}

    expression_type M_expr;
};
/**
   \ingroup vf

   \code
   auto e = Px();
   auto plus_e = +e;
   \endcode

   \return (unary) plus an expression
*/

template <class T> inline
Expr< UnaryPlus< Expr<T> > >
operator + ( const Expr<T>& expr )
{
    typedef UnaryPlus< Expr<T> > expr_t;

    return Expr< expr_t >( expr_t( expr ) );
}

/**
   \class UnaryMinus
   \brief handler for unary minus expression
*/
template < class ExprT >
class UnaryMinus : public ExprDynamicBase
{
public:

    static const size_type context = ExprT::context;
    static const bool is_terminal = false;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = ExprT::template HasTestFunction<Func>::result;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = ExprT::template HasTrialFunction<Func>::result;
    };
    template<typename Func>
    static const bool has_test_basis = ExprT::template has_test_basis<Func>;
    template<typename Func>
    static const bool has_trial_basis = ExprT::template has_trial_basis<Func>;
    using test_basis = typename ExprT::test_basis;
    using trial_basis = typename ExprT::trial_basis;

    typedef ExprT expression_type;
    typedef typename ExprT::value_type value_type;
    using evaluate_type = typename expression_type::evaluate_type;
    typedef UnaryMinus<ExprT> this_type;

    UnaryMinus( const ExprT& expr )
        :
        M_expr( expr )
    {
        ;
    }

    UnaryMinus( UnaryMinus const& u )
        :
        M_expr( u.M_expr )
    {
        ;
    }

    UnaryMinus( UnaryMinus&& u )
        :
        M_expr( u.M_expr )
    {
    }

    UnaryMinus&
    operator=( UnaryMinus const& u )
    {
        if ( this != &u )
            M_expr = u.M_expr;
        return *this;
    }

    template<typename... TheExpr>
    struct Lambda
    {
        typedef UnaryMinus<typename ExprT::template Lambda<TheExpr...>::type> type;
    };

    template<typename... TheExpr>
    typename Lambda<TheExpr...>::type
    operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_expr( e... ) ); }

    size_type dynamicContext() const { return vf::dynamicContext( M_expr ); } 

    //! polynomial order
    uint16_type polynomialOrder() const { return M_expr.polynomialOrder(); }

    //! expression is polynomial?
    bool isPolynomial() const { return M_expr.isPolynomial(); }

    //! evaluate the expression without context
    evaluate_type evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const
        {
            return -M_expr.evaluate(p,worldcomm);
        }

    expression_type const& expression() const
    {
        return M_expr;
    }

    void setParameterValues( std::map<std::string,value_type> const& mp )
        {
            M_expr.setParameterValues( mp );
        }
    void updateParameterValues( std::map<std::string,double> & pv ) const
        {
            M_expr.updateParameterValues( pv );
        }

    template <typename SymbolsExprType>
    auto applySymbolsExpr( SymbolsExprType const& se ) const
        {
            auto newExpr = M_expr.applySymbolsExpr( se );
            using new_expr_type = std::decay_t<decltype(newExpr)>;
            return UnaryMinus<new_expr_type>( newExpr );
        }

    template <typename TheSymbolExprType>
    bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const
        {
            return M_expr.hasSymbolDependency( symb, se );
        }
    template <typename TheSymbolExprType>
    void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const
        {
            M_expr.dependentSymbols( symb, res, se );
        }

    template <int diffOrder, typename TheSymbolExprType>
    auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr,
               TheSymbolExprType const& se ) const
        {
            auto theDiffExpr = M_expr.template diff<diffOrder>( diffVariable, world, dirLibExpr, se );
            using new_expr_type = std::decay_t<decltype(theDiffExpr)>;
            return UnaryMinus<new_expr_type>( theDiffExpr );
        }


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t>
    struct tensor
    {
        typedef typename expression_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_type;
        typedef typename tensor_expr_type::value_type value_type;
        typedef typename tensor_expr_type::shape shape;

        struct is_zero
        {
            static const bool value = tensor_expr_type::is_zero::value;
        };

        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev,
                Basis_j_t const& feu )
            :
            M_t_expr( expr.expression(), geom, fev, feu )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom,
                Basis_i_t const& fev )
            :
            M_t_expr( expr.expression(), geom, fev )
        {}
        tensor( this_type const& expr,
                Geo_t const& geom )
            :
            M_t_expr( expr.expression(), geom )
        {}
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs )
            :
            M_t_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... )
            {}
        template<typename IM>
        void init( IM const& im )
        {
            M_t_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev , Basis_j_t const& feu )
        {
            M_t_expr.update( geom, fev, feu );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev  )
        {
            M_t_expr.update( geom, fev );
        }
        void update( Geo_t const& geom )
        {
            M_t_expr.update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_t_expr.update( geom, face );
        }
        template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType>
        void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse,
                     Geo_t const& geom, const TheArgsType&... theUpdateArgs )
            {
                M_t_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs... );
            }

        template<typename ... CTX>
        void updateContext( CTX const& ... ctx )
        {
            M_t_expr.updateContext( ctx... );
        }

        value_type
        evalij( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2 ) const
        {
            return -M_t_expr.evalij( i, j, c1, c2 );
        }

        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -M_t_expr.evalijq( i, j, c1, c2, q );
        }
        template<int PatternContext>
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q,
                 mpl::int_<PatternContext> ) const
        {
            return -M_t_expr.evalijq( i, j, c1, c2, q, mpl::int_<PatternContext>() );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -M_t_expr.evaliq( i, c1, c2, q );
        }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return -M_t_expr.evalq( c1, c2, q );
        }
        tensor_expr_type M_t_expr;
    };
protected:
    UnaryMinus() {}

    expression_type M_expr;
};

/**
   \ingroup vf

   \code
   auto e = Px();
   auto minus_e = -e;
   \endcode

   \return (unary) minus an expression
*/
template <class T> inline
Expr< UnaryMinus< Expr<T> > >
operator - ( const Expr<T>& expr )
{
    typedef UnaryMinus< Expr<T> > expr_t;

    return Expr< expr_t >( expr_t( expr ) );
}
} // vf
} // Feel
#endif /* FEELPP_VF_UNARY_HPP */
