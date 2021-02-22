/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-06-07

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2012 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2019 Université de Strasbourg

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
   \file functions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-06-07
 */
#if !defined( STD_MATH_UNARY_FUNCTORS_HPP )
#define STD_MATH_UNARY_FUNCTORS_HPP 1

# include <boost/preprocessor/comparison/less.hpp>
# include <boost/preprocessor/logical/and.hpp>
# include <boost/preprocessor/control/if.hpp>
# include <boost/preprocessor/list/at.hpp>
# include <boost/preprocessor/list/cat.hpp>
# include <boost/preprocessor/list/for_each_product.hpp>
# include <boost/preprocessor/logical/or.hpp>
# include <boost/preprocessor/tuple/to_list.hpp>
# include <boost/preprocessor/tuple/eat.hpp>
# include <boost/preprocessor/facilities/empty.hpp>
# include <boost/preprocessor/punctuation/comma.hpp>
# include <boost/preprocessor/facilities/identity.hpp>

#include <boost/utility/enable_if.hpp>

/// \cond detail
#include <feel/feelcore/traits.hpp>
#include <feel/feelvf/unaryfunctor.hpp>
#include <feel/feelvf/binaryfunctor.hpp>
#include <feel/feelvf/arithmetic.hpp>

namespace Feel
{
namespace vf
{
namespace details
{
/**
   sign function
 */
template<typename T>
T
sign( T const& x )
{
    return ( x > T( 0.0 ) )?T( 1.0 ):( ( x < T( 0.0 ) )?T( -1.0 ):0.0 );
}

} // details

}
}

# /* Information about functions  */
#
# /* Accessors for the operator datatype. */
# define VF_FUNC_SYMBOL(O)        BOOST_PP_TUPLE_ELEM(8, 0, O)
# define VF_FUNC_NAME(O)          BOOST_PP_TUPLE_ELEM(8, 1, O)
# define VF_FUNC_IMPL(O)          BOOST_PP_TUPLE_ELEM(8, 2, O)
# define VF_FUNC_NAME_STRING(O)   BOOST_PP_TUPLE_ELEM(8, 3, O)
# define VF_FUNC_DOMAIN(O)        BOOST_PP_TUPLE_ELEM(8, 4, O)
# define VF_FUNC_IS_FLOATING(O )  BOOST_PP_TUPLE_ELEM(8, 5, O)
# define VF_FUNC_IS_LOGICAL(O)    BOOST_PP_TUPLE_ELEM(8, 6, O)
# define VF_FUNC_POLYNOMIALORDER_SCALING(O) BOOST_PP_TUPLE_ELEM(8, 7, O)
#
# /* List of applicative unary functions. */
# define VF_APPLICATIVE_UNARY_FUNCS \
   BOOST_PP_TUPLE_TO_LIST( \
      17, \
      ( \
         ( abs  , __Abs__ , Feel::math::abs     ,"absolute value"      , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( cos  , __Cos__ , Feel::math::cos     ,"cosine"              , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( sin  , __Sin__ , Feel::math::sin     ,"sine"                , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( tan  , __Tan__ , Feel::math::tan     ,"tangent"             , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( acos , __ACos__, Feel::math::acos    ,"inverse cosine"      , BoundedDomain<value_type>(-1.0,1.0), 1, 0, 2), \
         ( asin , __ASin__, Feel::math::asin    ,"inverse sine"        , BoundedDomain<value_type>(-1.0,1.0), 1, 0, 2), \
         ( atan , __ATan__, Feel::math::atan    ,"inverse tangent"     , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( cosh , __Cosh__, Feel::math::cosh    ,"hyperbolic cosine"   , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( sinh , __Sinh__, Feel::math::sinh    ,"hyperbolic sine"     , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( tanh , __Tanh__, Feel::math::tanh    ,"hyperbolic tangent"  , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( exp  , __Exp__ , Feel::math::exp     ,"exponential"         , UnboundedDomain<value_type>()      , 1, 0, 2), \
         ( log  , __Log__ , Feel::math::log     ,"logarithm"           , PositiveDomain<value_type>()       , 1, 0, 2), \
         ( sqrt , __Sqrt__, Feel::math::sqrt    ,"square root"         , PositiveDomain<value_type>()       , 1, 0, 2), \
         ( floor, __Floor__, std::floor         ,"floor"               , UnboundedDomain<value_type>()      , 1, 0, 1), \
         ( ceil , __Ceil__, std::ceil           ,"ceil"                , UnboundedDomain<value_type>()      , 1, 0, 1), \
         ( sign , __Sign__, details::sign       ,"sign"                , UnboundedDomain<value_type>()      , 1, 0, 1), \
         ( chi  , __Chi__ ,                     ,"chi"                 , UnboundedDomain<value_type>()      , 0, 1, 1) \
      ) \
   ) \
   /**/
# /* List of applicative unary functions. */
# define VF_APPLICATIVE_BINARY_FUNCS \
   BOOST_PP_TUPLE_TO_LIST( \
      1, \
      ( \
       ( atan2  , __ATan2__ , Feel::math::atan2 ,"arctan(y/x)" , std::make_pair(UnboundedDomain<value_type>(),UnboundedDomain<value_type>()) , 1, 0, 2) \
      ) \
   ) \
   /**/
#
#
# /* Generates code for all binary operators and integral type pairs. */
# define VF_UNARY_FUNCTIONS(_, O) \
      VF_UNARY_FUNCTIONS_CODE O \
   /**/
# define VF_BINARY_FUNCTIONS(_, O) \
      VF_BINARY_FUNCTIONS_CODE O \
   /**/



#define VF_UNARY_FUNCTIONS_CODE(O)                                      \
    template < typename ExprT1 >                                        \
    class VF_FUNC_NAME( O ) : public UnaryFunctor<typename ExprT1::value_type>, public ExprDynamicBase \
    {                                                                   \
    public:                                                             \
        using super2 = ExprDynamicBase;                                 \
        static const size_type context = ExprT1::context;               \
        static const bool is_terminal = false;                          \
                                                                        \
        template<typename Func>                                         \
            struct HasTestFunction                                      \
        {                                                               \
            static const bool result = false;                           \
        };                                                              \
                                                                        \
        template<typename Func>                                         \
            struct HasTrialFunction                                     \
        {                                                               \
            static const bool result = false;                           \
        };                                                              \
        template<typename Func>                                         \
            static const bool has_test_basis = false;                   \
        template<typename Func>                                         \
            static const bool has_trial_basis = false;                  \
        using test_basis = std::nullptr_t;                              \
        using trial_basis = std::nullptr_t;                             \
                                                                        \
        typedef UnaryFunctor<typename ExprT1::value_type> super;        \
        typedef typename super::functordomain_type functordomain_type;  \
        typedef typename super::functordomain_ptrtype functordomain_ptrtype; \
        typedef ExprT1 expression_1_type;                               \
        typedef VF_FUNC_NAME(O)<ExprT1> this_type;                      \
        typedef typename expression_1_type::value_type value_1_type;    \
        typedef value_1_type value_type;                                \
        using evaluate_type = typename expression_1_type::evaluate_type; \
                                                                        \
        VF_CHECK_ARITHMETIC_TYPE(value_1_type)                          \
                                                                        \
            explicit VF_FUNC_NAME(O)( expression_1_type const& __expr1  ) \
            :                                                           \
            super( VF_FUNC_NAME_STRING(O), functordomain_ptrtype(new VF_FUNC_DOMAIN(O) )), \
            super2( Feel::vf::dynamicContext( __expr1 ) ),              \
            M_expr_1( __expr1 )                                         \
            {                                                           \
                DVLOG(2) << "VF_FUNC_NAME(O)::VF_FUNC_NAME(O) default constructor\n"; \
            }                                                           \
                                                                        \
        VF_FUNC_NAME(O)( VF_FUNC_NAME(O) const& __vfp  )                \
            :                                                           \
            super( VF_FUNC_NAME_STRING(O), functordomain_ptrtype(new VF_FUNC_DOMAIN(O) )), \
            super2( Feel::vf::dynamicContext( __vfp.M_expr_1 ) ),       \
            M_expr_1( __vfp.M_expr_1 )                                \
                {                                                       \
                    DVLOG(2) << "VF_FUNC_NAME(O)::VF_FUNC_NAME(O) copy constructor\n"; \
                }                                                       \
                                                                        \
        bool isSymetric() const { return false; }                       \
                                                                        \
        uint16_type polynomialOrder() const { return VF_FUNC_POLYNOMIALORDER_SCALING(O)*M_expr_1.polynomialOrder(); } \
                                                                        \
        bool isPolynomial() const { return  ( M_expr_1.isPolynomial() && ( M_expr_1.polynomialOrder() == 0 ) ); } \
                                                                        \
        evaluate_type evaluate( bool p, worldcomm_ptr_t const& worldcomm ) const \
        {                                                               \
            auto eval = M_expr_1.evaluate(p,worldcomm);                 \
            evaluate_type res(eval.rows(),eval.cols());                 \
            for ( uint16_type i=0;i< eval.rows();++i )                  \
                for ( uint16_type j=0;j< eval.cols();++j )              \
                    res(i,j) = VF_FUNC_IMPL(O)( eval(i,j) );            \
            return res;                                                 \
        }                                                               \
                                                                        \
        void eval( int nx, value_type const* x, value_type* f ) const   \
        {                                                               \
            for( int i = 0; i < nx; ++i )                               \
                f[i] = VF_FUNC_IMPL(O)( x[i] );                         \
        }                                                               \
        template<typename... TheExpr>                                      \
        struct Lambda                                                   \
        {                                                               \
            typedef VF_FUNC_NAME( O )<typename ExprT1::template Lambda<TheExpr...>::type> type; \
        };                                                              \
                                                                        \
        template<typename... TheExpr>                                        \
            typename Lambda<TheExpr...>::type                           \
        operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_expr_1( e... ) ); } \
                                                                        \
                                                                        \
        expression_1_type const& expression() const { return M_expr_1; } \
                                                                        \
        void setParameterValues( std::map<std::string,value_type> const& mp ) \
        {                                                               \
            M_expr_1.setParameterValues( mp );                          \
        }                                                               \
        void updateParameterValues( std::map<std::string,double> & pv ) const \
        {                                                               \
            M_expr_1.updateParameterValues( pv );                       \
        }                                                               \
                                                                        \
        template <typename SymbolsExprType>                             \
            auto applySymbolsExpr( SymbolsExprType const& se ) const    \
        {                                                               \
            return VF_FUNC_SYMBOL( O ) ( M_expr_1.applySymbolsExpr( se ) ); \
        }                                                               \
                                                                        \
        template <typename TheSymbolExprType>                           \
            bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const \
        {                                                               \
            return M_expr_1.hasSymbolDependency( symb, se );            \
        }                                                               \
                                                                        \
        template <typename TheSymbolExprType>                           \
            void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const \
        {                                                               \
            M_expr_1.dependentSymbols( symb,res,se );                   \
        }                                                               \
                                                                        \
        template <int diffOrder, typename TheSymbolExprType>            \
            auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr, \
                       TheSymbolExprType const& se ) const              \
        {                                                               \
            CHECK( false ) << "TODO";                                   \
            return *this;                                               \
        }                                                               \
                                                                        \
        template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
            struct tensor                                               \
        {                                                               \
            typedef this_type expression_type;                          \
            typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t,Basis_j_t> tensor2_expr_type; \
            typedef typename tensor2_expr_type::value_type value_type;  \
            typedef typename vf::detail::ExtractGm<Geo_t>::gmc_ptrtype gmc_ptrtype; \
            typedef typename vf::detail::ExtractGm<Geo_t>::gmc_type gmc_type; \
            typedef typename tensor2_expr_type::shape shape;            \
                                                                        \
            struct is_zero { static const bool value = tensor2_expr_type::is_zero::value; }; \
                                                                        \
            tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                :                                                       \
                M_expr( expr.expression(), geom ),                     \
                M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        /*update( geom );*/                             \
                    }                                                   \
            tensor( this_type const& expr,Geo_t const& geom, Basis_i_t const& /*fev*/ ) \
                :                                                       \
                M_expr( expr.expression(), geom ),                     \
                M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        /*update( geom );*/                             \
                    }                                                   \
            tensor( this_type const& expr, Geo_t const& geom )             \
                :                                                       \
                M_expr( expr.expression(), geom ),                     \
                M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )         \
                    {                                                   \
                        /*update( geom );*/                             \
                    }                                                   \
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                        this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs ) \
                :                                                       \
                M_expr( std::true_type{}, exprExpanded.expression(), ttse, expr.expression(), geom, theInitArgs... ), \
                M_gmc( vf::detail::ExtractGm<Geo_t>::get( geom ) )      \
            {}                                                          \
                                                                        \
            template<typename IM>                                       \
                void init( IM const& im )                               \
            {                                                           \
                M_expr.init( im );                                     \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
            {                                                           \
                update( geom );                                         \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )  \
            {                                                           \
                update( geom );                                         \
            }                                                           \
            void update( Geo_t const& geom )                            \
            {                                                           \
                M_expr.update( geom );                                 \
            }                                                           \
            void update( Geo_t const& geom, uint16_type face )          \
            {                                                           \
                M_expr.update( geom, face );                           \
            }                                                           \
            template<typename ... CTX>                                  \
                void updateContext( CTX const& ... ctx )                \
            {                                                           \
                M_expr.updateContext( ctx... );                         \
            }                                                           \
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                             Geo_t const& geom, const TheArgsType&... theUpdateArgs ) \
            {                                                           \
                M_expr.update( std::true_type{}, exprExpanded.expression(), ttse, geom, theUpdateArgs...); \
            }                                                           \
                                                                        \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            template<int PatternContext> \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q, \
                         mpl::int_<PatternContext> ) const              \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
                                               \
                value_type                                              \
                evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return VF_FUNC_IMPL(O)( M_expr.evalq( c1, c2, q ) );   \
            }                                                           \
        private:                                                        \
            tensor2_expr_type M_expr;                                  \
            gmc_ptrtype M_gmc;                                         \
        };                                                              \
                                                                        \
    protected:                                                          \
        VF_FUNC_NAME(O)() {}                                            \
                                                                        \
        expression_1_type M_expr_1;                                    \
    };                                                                  \
                                                                        \
    template<typename ExprT1>                                           \
    inline                                                              \
    Expr< VF_FUNC_NAME( O )<typename mpl::if_<boost::is_arithmetic<ExprT1>, \
                                              mpl::identity<Cst<ExprT1> >, \
                                              mpl::identity<Expr<ExprT1> > >::type::type > > \
    VF_FUNC_SYMBOL( O )( Expr<ExprT1> const& __e1 )                     \
    {                                                                   \
        typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,         \
            mpl::identity<Cst<ExprT1> >,                                \
            mpl::identity<Expr<ExprT1> > >::type::type t1;               \
        typedef VF_FUNC_NAME(O)<t1> expr_t;                             \
        return Expr< expr_t >(  expr_t( t1( __e1 ) ) );                 \
    }                                                                   \
    template<typename ExprT1>                                           \
    inline                                                              \
    Expr< VF_FUNC_NAME( O )<Cst<ExprT1> > >                             \
    VF_FUNC_SYMBOL( O )( ExprT1 const& __e1, typename boost::enable_if<boost::is_arithmetic<ExprT1> >::type* dummy = 0 ) \
    {                                                                   \
        typedef Cst<ExprT1> t1;                                         \
        typedef VF_FUNC_NAME(O)<t1> expr_t;                             \
        return Expr< expr_t >(  expr_t( t1( __e1 ) ) );                 \
    }                                                                   \
    /**/


#define VF_BINARY_FUNCTIONS_CODE(O)                                     \
    template < typename ExprT1,typename ExprT2 >                                       \
    class VF_FUNC_NAME( O ) : public BinaryFunctor<typename ExprT1::value_type,typename ExprT2::value_type> \
    {                                                                   \
    public:                                                             \
                                                                        \
        static const size_type context = ExprT1::context|ExprT2::context; \
        static const bool is_terminal = false;                          \
                                                                        \
        template<typename Func>                                         \
            struct HasTestFunction                                      \
        {                                                               \
            static const bool result = false;                           \
        };                                                              \
                                                                        \
        template<typename Func>                                         \
            struct HasTrialFunction                                     \
        {                                                               \
            static const bool result = false;                           \
        };                                                              \
        using test_basis = std::nullptr_t;                              \
        using trial_basis = std::nullptr_t;                             \
                                                                        \
        typedef BinaryFunctor<typename ExprT1::value_type,typename ExprT2::value_type> super;      \
        typedef typename super::functordomain_1_type functordomain_1_type;  \
        typedef typename super::functordomain_2_type functordomain_2_type;  \
        /*typedef typename super::functordomain_ptrtype functordomain_ptrtype;*/ \
        typedef ExprT1 expression_1_type;                               \
        typedef ExprT2 expression_2_type;                               \
        typedef VF_FUNC_NAME(O)<ExprT1,ExprT2> this_type;               \
        typedef typename expression_1_type::value_type value_1_type;    \
        typedef typename expression_2_type::value_type value_2_type;    \
        typedef value_1_type value_type;                                \
        using evaluate_type = typename expression_1_type::evaluate_type; \
                                                                        \
        VF_CHECK_ARITHMETIC_TYPE(value_1_type)                          \
        VF_CHECK_ARITHMETIC_TYPE(value_2_type)                          \
                                                                        \
            explicit VF_FUNC_NAME(O)( expression_1_type const& __expr1, expression_2_type const& __expr2  ) \
            :                                                           \
            super( VF_FUNC_NAME_STRING(O), VF_FUNC_DOMAIN(O) ),         \
            M_expr_1( __expr1 ),                                        \
            M_expr_2( __expr2 )                                         \
            {                                                           \
                DVLOG(2) << "VF_FUNC_NAME(O)::VF_FUNC_NAME(O) default constructor\n"; \
            }                                                           \
                                                                        \
        VF_FUNC_NAME(O)( VF_FUNC_NAME(O) const& __vfp  )                \
            :                                                           \
            super( VF_FUNC_NAME_STRING(O), VF_FUNC_DOMAIN(O) ),         \
            M_expr_1( __vfp.M_expr_1 ),                                 \
            M_expr_2( __vfp.M_expr_2 )                                  \
                {                                                       \
                    DVLOG(2) << "VF_FUNC_NAME(O)::VF_FUNC_NAME(O) copy constructor\n"; \
                }                                                       \
                                                                        \
        bool isSymetric() const { return false; }                       \
                                                                        \
        uint16_type polynomialOrder() const { return VF_FUNC_POLYNOMIALORDER_SCALING(O)*std::max(M_expr_1.polynomialOrder(),M_expr_2.polynomialOrder()); } \
                                                                        \
        bool isPolynomial() const { return  ( M_expr_1.isPolynomial() && M_expr_2.isPolynomial() && ( M_expr_1.polynomialOrder() == 0 ) && ( M_expr_2.polynomialOrder() == 0 ) ); } \
                                                                        \
        evaluate_type evaluate( bool p, worldcomm_ptr_t const& worldcomm ) const \
        {                                                               \
            auto eval1 = M_expr_1.evaluate(p,worldcomm);                \
            auto eval2 = M_expr_2.evaluate(p,worldcomm);                \
            CHECK( eval1.rows() == eval2.rows() && eval1.cols() == eval2.cols() ) << "should be have same shape"; \
            evaluate_type res(eval1.rows(),eval1.cols());               \
            for ( uint16_type i=0;i< eval1.rows();++i )                 \
                for ( uint16_type j=0;j< eval1.cols();++j )             \
                    res(i,j) = VF_FUNC_IMPL(O)( eval1(i,j),eval2(i,j) ); \
            return res;                                                 \
        }                                                               \
                                                                        \
        void eval( int nx, value_1_type const* x, value_2_type const* y, value_type* f ) const   \
        {                                                               \
            for( int i = 0; i < nx; ++i )                               \
                f[i] = VF_FUNC_IMPL(O)( x[i],y[i] );                    \
        }                                                               \
        template<typename... TheExpr>                                   \
        struct Lambda                                                   \
        {                                                               \
            typedef VF_FUNC_NAME( O )<typename ExprT1::template Lambda<TheExpr...>::type, typename ExprT2::template Lambda<TheExpr...>::type > type; \
        };                                                              \
                                                                        \
        template<typename... TheExpr>                                        \
            typename Lambda<TheExpr...>::type                               \
        operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_expr_1( e... ), M_expr_2( e... ) ); } \
                                                                        \
                                                                        \
        expression_1_type const& expression1() const { return M_expr_1; } \
        expression_2_type const& expression2() const { return M_expr_2; } \
                                                                        \
        void setParameterValues( std::map<std::string,value_type> const& mp ) \
        {                                                               \
            M_expr_1.setParameterValues( mp );                          \
            M_expr_2.setParameterValues( mp );                          \
        }                                                               \
        void updateParameterValues( std::map<std::string,double> & pv ) const \
        {                                                               \
            M_expr_1.updateParameterValues( pv );                       \
            M_expr_2.updateParameterValues( pv );                       \
        }                                                               \
                                                                        \
        template <typename SymbolsExprType>                             \
            auto applySymbolsExpr( SymbolsExprType const& se ) const    \
        {                                                               \
            return VF_FUNC_SYMBOL( O ) ( M_expr_1.applySymbolsExpr( se ), M_expr_2.applySymbolsExpr( se ) ); \
        }                                                               \
                                                                        \
        template <typename TheSymbolExprType>                           \
            bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const \
        {                                                               \
            return M_expr_1.hasSymbolDependency( symb, se ) || M_expr_2.hasSymbolDependency( symb, se ); \
        }                                                               \
                                                                        \
        template <typename TheSymbolExprType>                           \
            void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const \
        {                                                               \
            M_expr_1.dependentSymbols( symb,res,se );                   \
            M_expr_2.dependentSymbols( symb,res,se );                   \
        }                                                               \
                                                                        \
        template <int diffOrder, typename TheSymbolExprType>            \
            auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr, \
                       TheSymbolExprType const& se ) const              \
        {                                                               \
            CHECK( false ) << "TODO";                                   \
            return *this;                                               \
        }                                                               \
                                                                        \
        template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
            struct tensor                                               \
        {                                                               \
            typedef this_type expression_type;                          \
            typedef typename expression_1_type::template tensor<Geo_t, Basis_i_t,Basis_j_t> tensor2_expr_1_type; \
            typedef typename expression_2_type::template tensor<Geo_t, Basis_i_t,Basis_j_t> tensor2_expr_2_type; \
            typedef typename tensor2_expr_1_type::value_type value_type;  \
            typedef typename vf::detail::ExtractGm<Geo_t>::gmc_ptrtype gmc_ptrtype; \
            typedef typename vf::detail::ExtractGm<Geo_t>::gmc_type gmc_type; \
            typedef typename tensor2_expr_1_type::shape shape;          \
                                                                        \
            struct is_zero { static const bool value = false/*tensor2_expr_type::is_zero::value*/; }; \
                                                                        \
            tensor( this_type const& expr, Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
                :                                                       \
                M_expr1( expr.expression1(), geom ),                    \
                M_expr2( expr.expression2(), geom ),                    \
                M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )       \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            tensor( this_type const& expr,Geo_t const& geom, Basis_i_t const& /*fev*/ ) \
                :                                                       \
                M_expr1( expr.expression1(), geom ),                    \
                M_expr2( expr.expression2(), geom ),                    \
                M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )       \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            tensor( this_type const& expr, Geo_t const& geom )          \
                :                                                       \
                M_expr1( expr.expression1(), geom ),                    \
                M_expr2( expr.expression2(), geom ),                    \
                M_gmc(vf::detail::ExtractGm<Geo_t>::get( geom ) )       \
                    {                                                   \
                        update( geom );                                 \
                    }                                                   \
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                        this_type const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs ) \
                :                                                       \
                M_expr1( std::true_type{}, exprExpanded.expression1(), ttse, expr.expression1(), geom, theInitArgs... ), \
                M_expr2( std::true_type{}, exprExpanded.expression2(), ttse, expr.expression2(), geom, theInitArgs... ), \
                M_gmc( vf::detail::ExtractGm<Geo_t>::get( geom ) )      \
                {}                                                      \
                                                                        \
            template<typename IM>                                       \
                void init( IM const& im )                               \
            {                                                           \
                M_expr1.init( im );                                     \
                M_expr2.init( im );                                     \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) \
            {                                                           \
                update( geom );                                         \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )  \
            {                                                           \
                update( geom );                                         \
            }                                                           \
            void update( Geo_t const& geom )                            \
            {                                                           \
                M_expr1.update( geom );                                 \
                M_expr2.update( geom );                                 \
            }                                                           \
            void update( Geo_t const& geom, uint16_type face )          \
            {                                                           \
                M_expr1.update( geom, face );                           \
                M_expr2.update( geom, face );                           \
            }                                                           \
            template<typename ... CTX>                                  \
                void updateContext( CTX const& ... ctx )                \
            {                                                           \
                M_expr1.updateContext( ctx... );                        \
                M_expr2.updateContext( ctx... );                        \
            }                                                           \
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                             Geo_t const& geom, const TheArgsType&... theUpdateArgs ) \
            {                                                           \
                M_expr1.update( std::true_type{}, exprExpanded.expression1(), ttse, geom, theUpdateArgs...); \
                M_expr2.update( std::true_type{}, exprExpanded.expression2(), ttse, geom, theUpdateArgs...); \
            }                                                           \
                                                                        \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            template<int PatternContext> \
                value_type                                              \
                evalijq( uint16_type /*i*/, uint16_type /*j*/, uint16_type c1, uint16_type c2, uint16_type q, \
                         mpl::int_<PatternContext> ) const              \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
                                                                        \
                value_type                                              \
                evaliq( uint16_type /*i*/, uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return evalq( c1, c2, q );                              \
            }                                                           \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q ) const \
            {                                                           \
                return VF_FUNC_IMPL(O)( M_expr1.evalq( c1, c2, q ), M_expr2.evalq( c1, c2, q ) ); \
            }                                                           \
        private:                                                        \
            tensor2_expr_1_type M_expr1;                                \
            tensor2_expr_2_type M_expr2;                                \
            gmc_ptrtype M_gmc;                                          \
        };                                                              \
                                                                        \
    protected:                                                          \
        VF_FUNC_NAME(O)() {}                                            \
                                                                        \
        expression_1_type M_expr_1;                                     \
        expression_2_type M_expr_2;                                     \
    };                                                                  \
                                                                        \
    template<typename ExprT1,typename ExprT2>                           \
    inline                                                              \
    Expr< VF_FUNC_NAME( O )<typename mpl::if_<boost::is_arithmetic<ExprT1>, \
                                              mpl::identity<Cst<ExprT1> >, \
                                              mpl::identity<Expr<ExprT1> > >::type::type, \
                            typename mpl::if_<boost::is_arithmetic<ExprT2>, \
                                              mpl::identity<Cst<ExprT2> >, \
                                              mpl::identity<Expr<ExprT2> > >::type::type > >                                                                     \
    VF_FUNC_SYMBOL( O )( Expr<ExprT1> const& __e1, Expr<ExprT2> const& __e2 ) \
    {                                                                   \
        typedef typename mpl::if_<boost::is_arithmetic<ExprT1>,         \
                                  mpl::identity<Cst<ExprT1> >,          \
                                  mpl::identity<Expr<ExprT1> > >::type::type t1; \
        typedef typename mpl::if_<boost::is_arithmetic<ExprT2>,         \
                                  mpl::identity<Cst<ExprT2> >,          \
                                  mpl::identity<Expr<ExprT2> > >::type::type t2; \
        typedef VF_FUNC_NAME(O)<t1,t2> expr_t;                          \
        return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );     \
    }                                                                   \
    template<typename ExprT1,typename ExprT2>                           \
    inline                                                              \
    Expr< VF_FUNC_NAME( O )<Cst<ExprT1>,Cst<ExprT2> > >                            \
    VF_FUNC_SYMBOL( O )( ExprT1 const& __e1, ExprT2 const& __e2, typename boost::enable_if<boost::is_arithmetic<ExprT1> >::type* dummy = 0 ) \
    {                                                                   \
        typedef Cst<ExprT1> t1;                                         \
        typedef Cst<ExprT2> t2;                                         \
        typedef VF_FUNC_NAME(O)<t1,t2> expr_t;                          \
        return Expr< expr_t >(  expr_t( t1( __e1 ), t2( __e2 ) ) );     \
    }                                                                   \
    /**/
#

namespace Feel
{
namespace vf
{

BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_UNARY_FUNCTIONS, 1, ( VF_APPLICATIVE_UNARY_FUNCS ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_FUNCTIONS, 1, ( VF_APPLICATIVE_BINARY_FUNCS ) )
}
}
/// \endcond
#endif /* STD_MATH_UNARY_FUNCTORS_HPP */
