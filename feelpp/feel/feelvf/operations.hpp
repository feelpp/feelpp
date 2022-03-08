/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-06-01

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2011 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2011-2014 Feel++ Consortium

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
   \file operations.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-06-01
 */
#if !defined( FEELPP_OPERATIONS_HPP )
#define FEELPP_OPERATIONS_HPP 1

#include <feel/feelconfig.h>
#if defined( FEELPP_HAS_QD_H )
# include <qd/qd.h>
#endif
/// \cond detail
# include <complex>
# include <boost/preprocessor/comparison/less.hpp>
# include <boost/preprocessor/comparison/equal.hpp>
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
# include <feel/feelvf/pow.hpp>

# /* Information about C operators */
#
# /* Accessors for the operator datatype. */
# define VF_OP_SYMBOL(O)      BOOST_PP_TUPLE_ELEM(9, 0, O)
# define VF_OP_NAME(O)        BOOST_PP_TUPLE_ELEM(9, 1, O)
# define VF_OP_IS_FLOATING(O) BOOST_PP_TUPLE_ELEM(9, 2, O)
# define VF_OP_IS_LOGICAL(O)  BOOST_PP_TUPLE_ELEM(9, 3, O)
# define VF_OP_IS_SHIFT(O)    BOOST_PP_TUPLE_ELEM(9, 4, O)
# define VF_OP_IS_ADD(O)      BOOST_PP_TUPLE_ELEM(9, 5, O)
# define VF_OP_IS_SYMETRIC(O) BOOST_PP_TUPLE_ELEM(9, 6, O)
# define VF_OP_SHAPE(O)       BOOST_PP_TUPLE_ELEM(9, 7, O)
# define VF_OP_IMORDER(O)     BOOST_PP_TUPLE_ELEM(9, 8, O)
#
# /* List of applicative unary operators. */
# define VF_APPLICATIVE_UNARY_OPS \
   BOOST_PP_TUPLE_TO_LIST( \
      3, \
      ( \
         ( ! , vf_logical_not, 1, 1, 0,0,0), \
         ( ~ , vf_bitwise_not, 0, 0, 0,0,0), \
         ( - , vf_neg,         1, 0, 0,0,0) \
      ) \
   ) \
   /**/
#
#if 0
# /* List of applicative binary operators. */
# define VF_APPLICATIVE_BINARY_OPS \
   BOOST_PP_TUPLE_TO_LIST( \
      18, \
      ( \
         ( *  , vf_mul           ,1 ,0 ,0 ,0 ,1), \
         ( /  , vf_div           ,1 ,0 ,0 ,0 ,0), \
         ( %  , vf_mod           ,0 ,0 ,0 ,0 ,0), \
         ( +  , vf_add           ,1 ,0 ,0 ,1 ,1), \
         ( -  , vf_sub           ,1 ,0 ,0 ,1 ,0), \
         ( << , vf_shift_left    ,0 ,0 ,1 ,0 ,0), \
         ( >> , vf_shift_right   ,0 ,0 ,1 ,0 ,0), \
         ( <  , vf_less          ,1 ,1 ,0 ,0 ,0), \
         ( <= , vf_less_equal    ,1 ,1 ,0 ,0 ,0), \
         ( >= , vf_greater_equal ,1 ,1 ,0 ,0 ,0), \
         ( >  , vf_greater       ,1 ,1 ,0 ,0 ,0), \
         ( == , vf_equal         ,1 ,1 ,0 ,0 ,0), \
         ( != , vf_not_equal     ,1 ,1 ,0 ,0 ,0), \
         ( &  , vf_bitwise_and   ,0 ,0 ,0 ,0 ,0), \
         ( |  , vf_bitwise_or    ,0 ,0 ,0 ,0 ,0), \
         ( ^  , vf_bitwise_xor   ,0 ,0 ,0 ,0 ,0), \
         ( && , vf_logical_and   ,1 ,1 ,0 ,0 ,0), \
         ( || , vf_logical_or    ,1 ,1 ,0 ,0, 0) \
      ) \
   ) \
   /**/
#else
# // get rid of operator^ which is used for power computation
# /* List of applicative binary operators. */
# define VF_APPLICATIVE_BINARY_OPS \
   BOOST_PP_TUPLE_TO_LIST( \
      13, \
      ( \
         ( *  , vf_mul           ,1 ,0 ,0 ,0 ,1, shape_op_mul      ,2 ), \
         ( /  , vf_div           ,1 ,0 ,0 ,0 ,0, shape_op_div      ,1 ), \
         ( +  , vf_add           ,1 ,0 ,0 ,1 ,1, shape_op_samerank ,1 ), \
         ( %  , vf_mod           ,1 ,0 ,0 ,0 ,0, shape_op_div      ,1 ), \
         ( -  , vf_sub           ,1 ,0 ,0 ,1 ,0, shape_op_samerank ,1 ), \
         ( <  , vf_less          ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( <= , vf_less_equal    ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( >= , vf_greater_equal ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( >  , vf_greater       ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( == , vf_equal         ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( != , vf_not_equal     ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( && , vf_logical_and   ,1 ,1 ,0 ,0 ,0, shape_op_id       ,0 ), \
         ( || , vf_logical_or    ,1 ,1 ,0 ,0, 0, shape_op_id       ,0 )  \
      ) \
   ) \
   /**/
#endif
#
# /* Accessors for the type datatype. */
# define VF_TYPE_NAME(T)         BOOST_PP_TUPLE_ELEM(5, 0, T)
# define VF_TYPE_ABBREVIATION(T) BOOST_PP_TUPLE_ELEM(5, 1, T)
# define VF_TYPE_IS_FLOATING(T)  BOOST_PP_TUPLE_ELEM(5, 2, T)
# define VF_TYPE_RANK(T)         BOOST_PP_TUPLE_ELEM(5, 3, T)
# define VF_TYPE_IS_EXPR(T)      BOOST_PP_TUPLE_ELEM(5, 4, T)
#
#if 0
# define VF_BUILTIN_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      12, \
      ( \
         ( signed char    ,sc, 0, 1,0), \
         ( char           ,ch, 0, 1,0), \
         ( unsigned char  ,uc, 0, 1,0), \
         ( short          ,ss, 0, 2,0), \
         ( unsigned short ,us, 0, 2,0), \
         VF_TYPE_INT, \
         ( unsigned       ,ui, 0, 4,0), \
         ( long           ,sl, 0, 5,0), \
         ( unsigned long  ,ul, 0, 6,0), \
         ( float          ,fl, 1, 7,0), \
         ( double         ,db, 1, 8,0), \
         ( long double    ,ld, 1, 9,0)  \
      ) \
   ) \
   /**/
#else
# if defined( DISABLE_FEELPP_HAS_QD_H ) && defined( DISABLE_FEELPP_HAS_MPFR )
#  define VF_BUILTIN_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      12, \
      ( \
         VF_TYPE_INT, \
         ( unsigned       ,ui, 0, 4,0), \
         ( long           ,sl, 0, 5,0), \
         ( unsigned long  ,ul, 0, 6,0), \
         ( float          ,fl, 1, 7,0), \
         ( double         ,db, 1, 8,0), \
         ( long double    ,ld, 1, 9,0), \
         ( std::complex<float> ,cfl, 1, 10,0), \
         ( std::complex<double> ,cld, 1, 11,0), \
         ( mp_type        ,mp, 1,12,0), \
         ( dd_real        ,dd, 1,13,0), \
         ( qd_real        ,dd, 1,14,0) \
      ) \
   ) \
   /**/
# elif defined (DISABLE_FEELPP_HAS_QD_H )
#  define VF_BUILTIN_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      11, \
      ( \
         VF_TYPE_INT, \
         ( unsigned       ,ui, 0, 4,0), \
         ( long           ,sl, 0, 5,0), \
         ( unsigned long  ,ul, 0, 6,0), \
         ( float          ,fl, 1, 7,0), \
         ( double         ,db, 1, 8,0), \
         ( long double    ,ld, 1, 9,0), \
         ( std::complex<float> ,cfl, 1, 10,0), \
         ( std::complex<double> ,cld, 1, 11,0), \
         ( dd_real        ,dd, 1,12,0), \
         ( qd_real        ,dd, 1,13,0) \
      ) \
   ) \
   /**/
# elif defined (DISABLE_FEELPP_HAS_MPFR )
#  define VF_BUILTIN_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      10, \
      ( \
         VF_TYPE_INT, \
         ( unsigned       ,ui, 0, 4,0), \
         ( long           ,sl, 0, 5,0), \
         ( unsigned long  ,ul, 0, 6,0), \
         ( float          ,fl, 1, 7,0), \
         ( double         ,db, 1, 8,0), \
         ( long double    ,ld, 1, 9,0), \
         ( std::complex<float> ,cfl, 1, 10,0), \
         ( std::complex<double> ,cld, 1, 11,0), \
         ( mp_type        ,mp, 1,12,0)  \
      ) \
   ) \
   /**/
#else
#  define VF_BUILTIN_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      9, \
      ( \
         VF_TYPE_INT, \
         ( unsigned       ,ui, 0, 4,0), \
         ( long           ,sl, 0, 5,0), \
         ( unsigned long  ,ul, 0, 6,0), \
         ( float          ,fl, 1, 7,0), \
         ( double         ,db, 1, 8,0), \
         ( long double    ,ld, 1, 9,0), \
         ( std::complex<float> ,cfl, 1, 10,0), \
         ( std::complex<double> ,cld, 1, 11,0) \
        )                                      \
   ) \
   /**/
# endif
#endif
#
# /* Type int is needed in some type computations. */
# define VF_TYPE_INT (int, si, 0, 3, 0)
#
# define VF_VALUE_TYPE(L)  \
  BOOST_PP_LIST_CAT(BOOST_PP_TUPLE_TO_LIST(2,(value_type_,VF_TYPE_ABBREVIATION(L)))) \
   /**/
#
# define VF_TYPE_CV(L)                                                                         \
   BOOST_PP_IF(VF_TYPE_IS_EXPR(L),                                                             \
               BOOST_PP_IDENTITY(const&),                                                   \
               BOOST_PP_EMPTY )()                                                           \
   /**/
#
# define VF_TYPE_VALUE_TYPE(L)                                                                 \
   BOOST_PP_IF(VF_TYPE_IS_EXPR(L),                                                             \
               typename VF_TYPE_NAME(L)::value_type,                                           \
               VF_TYPE_NAME( L ) )                                                             \
   /**/
#
# define VF_TYPE_TYPE(L)                                                                       \
   BOOST_PP_IF(VF_TYPE_IS_EXPR(L),                                                             \
               Expr<VF_TYPE_NAME(L)>,                                   \
               Expr<Cst<VF_TYPE_NAME( L )>> )                           \
   /**/
#
# define VF_TYPE_TYPE_EXPR(L)                                                                  \
   BOOST_PP_IF(VF_TYPE_IS_EXPR(L),                                                             \
               Expr<VF_TYPE_NAME(L)>,                                                          \
               VF_TYPE_NAME( L ) )                                                             \
   /**/
#
# define VF_TYPE_TYPE_EXPR_CST(L)                                                              \
   BOOST_PP_IF(VF_TYPE_IS_EXPR(L),                                                             \
               Expr<VF_TYPE_NAME(L)>,                                                          \
               VF_TYPE_NAME( L )/*Expr<Cst<VF_TYPE_NAME( L )>>*/ )      \
   /**/
#
# define VF_TYPE_TYPE_CST(L)                                                                  \
   BOOST_PP_IF(BOOST_PP_NOT(VF_TYPE_IS_EXPR(L)),                                              \
               BOOST_PP_IDENTITY(cst/*Cst<VF_TYPE_NAME( L )>*/),        \
               BOOST_PP_EMPTY )()                                                             \
   /**/
#
# define VF_TEMPLATE_TYPE(L)                       \
   BOOST_PP_IF(VF_TYPE_IS_EXPR(L),                 \
               class VF_TYPE_NAME( L ),            \
               L )                              \
    /**/
#
#define VF_TEXT(z, n, text) text
#
# define VF_SPECIALIZATION_IF_BUILTIN(L,R)                         \
  BOOST_PP_IF(BOOST_PP_OR( BOOST_PP_NOT( VF_TYPE_IS_EXPR(L) ),     \
                           BOOST_PP_NOT( VF_TYPE_IS_EXPR(R) ) ),   \
              BOOST_PP_IDENTITY(</*VF_TYPE_TYPE*/VF_TYPE_NAME(L)), \
              BOOST_PP_EMPTY)()                                    \
  BOOST_PP_IF(BOOST_PP_OR( BOOST_PP_NOT( VF_TYPE_IS_EXPR(L) ),     \
                           BOOST_PP_NOT( VF_TYPE_IS_EXPR(R) ) ),   \
              BOOST_PP_COMMA,                                      \
              BOOST_PP_EMPTY)()                                    \
  BOOST_PP_IF(BOOST_PP_OR( BOOST_PP_NOT( VF_TYPE_IS_EXPR(L) ),     \
                           BOOST_PP_NOT( VF_TYPE_IS_EXPR(R) ) ),   \
              BOOST_PP_IDENTITY(/*VF_TYPE_TYPE*/VF_TYPE_NAME(R)>), \
              BOOST_PP_EMPTY)()                                    \
    /**/

#
#
# /* List of expression types. */
# define VF_EXPRL_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      1, \
      ( \
         ( ExprL  ,L,   1, 10, 1) \
      ) \
   ) \
   /**/
# /* List of expression types. */
# define VF_EXPRR_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      1, \
      ( \
         ( ExprR  ,R,   1, 10, 1) \
      ) \
   ) \
   /**/
#
# define VF_BOOST_PP_EMPTY() \
   /**/


# /* Generates code for all binary operators and integral type pairs. */
# define VF_BINARY_ARRAY_OP_DECLARATION(_, OLR) \
      VF_BINARY_ARRAY_OP_CODE_DECLARATION OLR \
   /**/

#define VF_BINARY_ARRAY_OP_CODE_DECLARATION(O,L,R)                      \
    template <BOOST_PP_IF( VF_TYPE_IS_EXPR( L ),                        \
                           BOOST_PP_IDENTITY(class VF_TYPE_NAME(L)),    \
                           BOOST_PP_EMPTY                               \
                           )()                                          \
              BOOST_PP_IF( BOOST_PP_AND( VF_TYPE_IS_EXPR(L),            \
                                         VF_TYPE_IS_EXPR(R) ),          \
                           BOOST_PP_COMMA,                              \
                           BOOST_PP_EMPTY )()                           \
              BOOST_PP_IF( VF_TYPE_IS_EXPR( R ),                        \
                           BOOST_PP_IDENTITY(class VF_TYPE_NAME(R)),    \
                           BOOST_PP_EMPTY                               \
                           )()>                                         \
    class VF_OP_NAME( O )  VF_SPECIALIZATION_IF_BUILTIN( L, R ) ;       \
    /**/




#
# /* Generates code for all binary operators and integral type pairs. */
# define VF_BINARY_ARRAY_OP(_, OLR) \
      VF_BINARY_ARRAY_OP_CODE OLR \
   /**/

#define VF_BINARY_ARRAY_OP_CODE(O,L,R)                                  \
    template <BOOST_PP_IF( VF_TYPE_IS_EXPR( L ),                        \
                           BOOST_PP_IDENTITY(class VF_TYPE_NAME(L)),    \
                           BOOST_PP_EMPTY                               \
                           )()                                          \
              BOOST_PP_IF( BOOST_PP_AND( VF_TYPE_IS_EXPR(L),            \
                                         VF_TYPE_IS_EXPR(R) ),          \
                           BOOST_PP_COMMA,                              \
                           BOOST_PP_EMPTY )()                           \
              BOOST_PP_IF( VF_TYPE_IS_EXPR( R ),                        \
                           BOOST_PP_IDENTITY(class VF_TYPE_NAME(R)),    \
                           BOOST_PP_EMPTY                               \
                           )()>                                         \
    class VF_OP_NAME( O ) VF_SPECIALIZATION_IF_BUILTIN( L, R ) : public ExprDynamicBase        \
    {                                                                   \
      public:                                                           \
        typedef VF_OP_NAME( O )<VF_TYPE_NAME( L ), VF_TYPE_NAME( R )> expression_type; \
        using this_type = expression_type;                              \
        typedef VF_TYPE_TYPE( L ) L_type;                               \
        typedef VF_TYPE_TYPE( R ) R_type;                               \
        using value_left_type = typename L_type::value_type;            \
        using value_right_type = typename R_type::value_type;           \
                                                                        \
        static const bool is_op_mul = std::is_same_v< expression_type, vf_mul<VF_TYPE_NAME( L ), VF_TYPE_NAME( R )> >; \
        static const bool is_op_div = std::is_same_v< expression_type, vf_div<VF_TYPE_NAME( L ), VF_TYPE_NAME( R )> >; \
        static const bool is_op_add = std::is_same_v< expression_type, vf_add<VF_TYPE_NAME( L ), VF_TYPE_NAME( R )> >; \
        static const bool is_op_sub = std::is_same_v< expression_type, vf_sub<VF_TYPE_NAME( L ), VF_TYPE_NAME( R )> >; \
                                                                        \
        static const size_type context = L_type::context | R_type::context; \
        size_type dynamicContext() const { return vf::dynamicContext( M_left ) | vf::dynamicContext( M_right ); } \
                                                                        \
        static const bool is_terminal = false;                           \
                                                                        \
        template<typename Func>                                         \
            struct HasTestFunction                                      \
        {                                                               \
            static const bool result = L_type::template HasTestFunction<Func>::result|R_type::template HasTestFunction<Func>::result; \
        };                                                              \
                                                                        \
        template<typename Func>                                         \
            struct HasTrialFunction                                     \
        {                                                               \
            static const bool result = L_type::template HasTrialFunction<Func>::result|R_type::template HasTrialFunction<Func>::result; \
        };                                                              \
        template<typename Func>                                         \
            static const bool has_test_basis = L_type::template HasTestFunction<Func>::result|R_type::template HasTestFunction<Func>::result; \
        template<typename Func>                                         \
            static const bool has_trial_basis = L_type::template HasTrialFunction<Func>::result|R_type::template HasTrialFunction<Func>::result; \
        using test_basis =  typename mpl::if_< std::is_null_pointer<typename L_type::test_basis >, typename R_type::test_basis , typename L_type::test_basis >::type; \
        using trial_basis = std::nullptr_t;                             \
                                                                        \
        typedef typename mpl::if_<mpl::greater<mpl::sizeof_<value_left_type>, \
            mpl::sizeof_<value_right_type> >,                           \
                                  mpl::identity<value_left_type>,      \
                                  mpl::identity<value_right_type> >::type::type value_type; \
        using evaluate_type = Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic >; \
                                                                        \
        VF_OP_NAME( O )( L_type const& left, R_type const& right )      \
            :                                                           \
            M_left(left),                                              \
            M_right(right)                                             \
            {}                                                          \
        VF_OP_NAME( O )( VF_OP_NAME(O) const& __m )                     \
            :                                                           \
            M_left( __m.M_left ),                                     \
            M_right( __m.M_right )                                    \
                {;}                                                     \
        ~VF_OP_NAME( O )()                                              \
        {}                                                              \
        template<typename... TheExpr>                                   \
        struct Lambda                                                   \
        {                                                               \
            typedef VF_OP_NAME( O )<typename L_type::template Lambda<TheExpr...>::type::expression_type,typename R_type::template Lambda<TheExpr...>::type::expression_type> type; \
        };                                                              \
                                                                        \
        template<typename... TheExpr>                                      \
            typename Lambda<TheExpr...>::type                              \
            operator()( TheExpr... e ) { return typename Lambda<TheExpr...>::type( M_left(e...), M_right(e...) ); } \
                                                                        \
        uint16_type polynomialOrder() const {                           \
            uint8_type t = VF_OP_IMORDER(O);                            \
            if ( t == 0 ) return 0;                                     \
            else if ( t == 1 ) return std::max( M_left.polynomialOrder(), M_right.polynomialOrder() ); \
            else return M_left.polynomialOrder() + M_right.polynomialOrder(); \
        }                                                               \
        bool isPolynomial() const { return M_left.isPolynomial() && M_right.isPolynomial(); } \
                                                                       \
        L_type VF_TYPE_CV(L) left() const { return M_left; }           \
        R_type VF_TYPE_CV(R) right() const { return M_right; }         \
                                                                        \
        void setParameterValues( std::map<std::string,double> const& mp ) \
        {                                                               \
            M_left.setParameterValues( mp );                            \
            M_right.setParameterValues( mp );                           \
        }                                                               \
        void updateParameterValues( std::map<std::string,double> & pv ) const \
        {                                                               \
            M_left.updateParameterValues( pv );                         \
            M_right.updateParameterValues( pv );                        \
        }                                                               \
                                                                        \
        template <typename SymbolsExprType>                             \
            auto applySymbolsExpr( SymbolsExprType const& se ) const    \
        {                                                               \
            return M_left.applySymbolsExpr( se ) VF_OP_SYMBOL( O ) M_right.applySymbolsExpr( se ); \
        }                                                               \
                                                                        \
        template <typename TheSymbolExprType>                           \
            bool hasSymbolDependency( std::string const& symb, TheSymbolExprType const& se ) const \
        {                                                               \
            return M_left.hasSymbolDependency( symb, se ) || M_right.hasSymbolDependency( symb, se ); \
        }                                                               \
                                                                        \
        template <typename TheSymbolExprType>                           \
            void dependentSymbols( std::string const& symb, std::map<std::string,std::set<std::string>> & res, TheSymbolExprType const& se ) const \
        {                                                               \
            M_left.dependentSymbols( symb,res,se );                     \
            M_right.dependentSymbols( symb,res,se );                    \
        }                                                               \
                                                                        \
        template <int diffOrder, typename TheSymbolExprType>            \
            auto diff( std::string const& diffVariable, WorldComm const& world, std::string const& dirLibExpr, \
                       TheSymbolExprType const& se ) const              \
        {                                                               \
            auto ldiff = M_left.template diff<diffOrder>( diffVariable, world, dirLibExpr, se ); \
            auto rdiff = M_right.template diff<diffOrder>( diffVariable, world, dirLibExpr, se ); \
            if constexpr( is_op_mul )                                   \
                            return ldiff*M_right + M_left*rdiff;        \
            else if constexpr( is_op_div )                              \
                                 return (ldiff*M_right - M_left*rdiff)/pow(M_right,2); \
            else if constexpr( is_op_add || is_op_sub )                 \
                                 return ldiff VF_OP_SYMBOL( O ) rdiff;  \
            else { CHECK(false ) << "TODO or not possible";return *this; } \
        }                                                               \
                                                                        \
        template<typename Geo_t, typename Basis_i_t, typename Basis_j_t = Basis_i_t> \
            struct tensor                                               \
        {                                                               \
            typedef this_type expression_type;                          \
            typedef typename L_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> l_type; \
            typedef typename R_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> r_type; \
            typedef strongest_numeric_type<typename l_type::value_type, \
                                           typename r_type::value_type> value_type; \
            typedef typename VF_OP_SHAPE( O )<typename l_type::shape, typename r_type::shape>::type shape; \
            static const int shape_op = VF_OP_SHAPE( O )<typename l_type::shape, typename r_type::shape>::op; \
                                                                        \
            struct is_zero {                                            \
                static const bool value = VF_OP_SHAPE( O )<typename l_type::shape, typename r_type::shape>::template is_zero<l_type::is_zero::value, r_type::is_zero::value>::value; \
                static const bool update_and_eval_left = VF_OP_SHAPE( O )<typename l_type::shape, typename r_type::shape>::template is_zero<l_type::is_zero::value, r_type::is_zero::value>::update_and_eval_left; \
                static const bool update_and_eval_right = VF_OP_SHAPE( O )<typename l_type::shape, typename r_type::shape>::template is_zero<l_type::is_zero::value, r_type::is_zero::value>::update_and_eval_right; \
            };                                                          \
                                                                        \
            template<typename ExprT>                                    \
                tensor( ExprT const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) \
                :                                                       \
                M_left( expr.left(), geom, fev, feu ),                 \
                M_right( expr.right(), geom, fev, feu )                \
                    {                                                  \
                        DVLOG(2) << "Operation " BOOST_PP_STRINGIZE( VF_OP_SYMBOL( O ) ) " is_zero " << is_zero::value << " " \
                                      << "update_and_eval_left " << is_zero::update_and_eval_left << " " \
                                      << " update_and_eval_right " << is_zero::update_and_eval_right << " " \
                                      << " left_is_zero " << l_type::is_zero::value << " " \
                                      << " right_is_zero " << r_type::is_zero::value << "\n"; \
                    }                                                   \
            template<typename ExprT>                                    \
            tensor( ExprT const& expr,Geo_t const& geom, Basis_i_t const& fev ) \
                :                                                       \
                M_left( expr.left(), geom, fev ),                      \
                M_right( expr.right(), geom, fev )                     \
                    {}                                                  \
            template<typename ExprT>                                    \
                tensor( ExprT const& expr, Geo_t const& geom )          \
                :                                                       \
                M_left( expr.left(), geom ),                           \
                M_right( expr.right(), geom )                          \
                    {                                                   \
                    }                                                   \
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename ExprT, typename... TheArgsType> \
                tensor( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                        ExprT const& expr, Geo_t const& geom, const TheArgsType&... theInitArgs ) \
                :                                                       \
                M_left( std::true_type{}, exprExpanded.left(), ttse, expr.left(), geom, theInitArgs... ), \
                M_right( std::true_type{}, exprExpanded.right(), ttse, expr.right(), geom, theInitArgs... ) \
                {}                                                      \
                                                                        \
            void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) noexcept \
            {                                                           \
                if ( is_zero::update_and_eval_left )                    \
                    M_left.update( geom, fev, feu );                   \
                if ( is_zero::update_and_eval_right )                   \
                    M_right.update( geom, fev, feu );                  \
            }                                                           \
            void update( Geo_t const& geom, Basis_i_t const& fev ) noexcept \
            {                                                           \
                if ( is_zero::update_and_eval_left )                    \
                    M_left.update( geom, fev );                        \
                if ( is_zero::update_and_eval_right )                   \
                    M_right.update( geom, fev );                       \
            }                                                           \
            void update( Geo_t const& geom ) noexcept                   \
            {                                                           \
                if ( is_zero::update_and_eval_left )                    \
                    M_left.update( geom );                             \
                if ( is_zero::update_and_eval_right )                   \
                    M_right.update( geom );                            \
            }                                                           \
            template<typename ... CTX>                                  \
                void updateContext( CTX const& ... ctx ) noexcept       \
            {                                                           \
                M_left.updateContext( ctx... );                        \
                M_right.updateContext( ctx... );                        \
            }                                                           \
            template<typename TheExprExpandedType,typename TupleTensorSymbolsExprType, typename... TheArgsType> \
                void update( std::true_type /**/, TheExprExpandedType const& exprExpanded, TupleTensorSymbolsExprType & ttse, \
                             Geo_t const& geom, const TheArgsType&... theUpdateArgs ) \
            {                                                           \
                M_left.update( std::true_type{}, exprExpanded.left(), ttse, geom, theUpdateArgs... ); \
                M_right.update( std::true_type{}, exprExpanded.right(), ttse, geom, theUpdateArgs... ); \
            }                                                           \
                                                                        \
            value_type                                                  \
                evalij( uint16_type i, uint16_type j ) const noexcept   \
            {                                                           \
                return M_left.evalij(i,j) VF_OP_SYMBOL( O ) M_right.evalij(i,j); \
            }                                                           \
                              \
            value_type                                                  \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const noexcept\
            {                                                           \
                return evalijq(i,j,c1,c2,q,mpl::int_<shape_op>() );     \
            }                                                           \
            template<int PatternContext> \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext> ) const noexcept \
            {                                                           \
                return evalijq(i,j,c1,c2,q,mpl::int_<PatternContext>(),mpl::int_<shape_op>() ); \
            }                                                           \
                              \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0>  ) const noexcept \
            {                                                           \
                return evalijq( i, j, c1, c2, q, mpl::bool_<is_zero::update_and_eval_left>(), mpl::bool_<is_zero::update_and_eval_right>() ); \
            }                                                           \
                              \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<false>, mpl::bool_<false>  ) const noexcept \
            {                                                           \
                return value_type( 0 );                                 \
            }                                                           \
                              \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true>, mpl::bool_<false>  ) const noexcept \
            {                                                           \
                return M_left.evalijq(i, j, c1, c2, q);                \
            }                                                           \
                              \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<false>, mpl::bool_<true>  ) const noexcept \
            {                                                           \
                return M_right.evalijq(i, j, c1, c2, q);               \
            }                                                           \
                              \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::bool_<true>, mpl::bool_<true>  ) const noexcept \
            {                                                           \
                if constexpr ( l_type::shape::is_scalar && !r_type::shape::is_scalar ) \
                     return M_left.evalijq(i, j, 0, 0, q) VF_OP_SYMBOL( O ) M_right.evalijq(i,j, c1, c2, q); \
                else if constexpr ( !l_type::shape::is_scalar && r_type::shape::is_scalar ) \
                    return M_left.evalijq(i, j, c1, c2, q) VF_OP_SYMBOL( O ) M_right.evalijq(i,j, 0, 0, q); \
                else                                                    \
                    return M_left.evalijq(i, j, c1, c2, q) VF_OP_SYMBOL( O ) M_right.evalijq(i,j, c1, c2, q); \
            }                                                           \
                                                                        \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1>  ) const noexcept \
            {                                                           \
                if constexpr ( is_zero::value )                         \
                                 return value_type( 0 );                \
                else                                                    \
                {                                                       \
                    value_type res( value_type( 0 ) );                  \
                    for(uint16_type ii = 0; ii < l_type::shape::N; ++ii ) \
                        res += M_left.evalijq(i, j, c1, ii, q ) VF_OP_SYMBOL( O ) M_right.evalijq(i, j, ii, c2, q ); \
                    return res;                                         \
                }                                                       \
            }                                                           \
                                                                        \
            template<int PatternContext> \
                value_type                                              \
                evalijq__( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext>, mpl::int_<0>  ) const noexcept \
            {                                                           \
                if ( is_zero::value )                                   \
                    return value_type( 0 );                             \
                else if ( !is_zero::update_and_eval_left && is_zero::update_and_eval_right ) \
                    return M_right.evalijq(i,j, c1, c2, q,mpl::int_<PatternContext>()); \
                else if ( is_zero::update_and_eval_left && !is_zero::update_and_eval_right ) \
                    return M_left.evalijq(i, j, c1, c2, q,mpl::int_<PatternContext>()); \
                else                                                    \
                    return M_left.evalijq(i, j, c1, c2, q,mpl::int_<PatternContext>()) VF_OP_SYMBOL( O ) M_right.evalijq(i,j, c1, c2, q,mpl::int_<PatternContext>()); \
            }                                                           \
            template<int PatternContext> \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext>, mpl::int_<0>  ) const noexcept \
            {                                                           \
                return M_left.evalijq(i, j, c1, c2, q,mpl::int_<PatternContext>()) VF_OP_SYMBOL( O ) M_right.evalijq(i,j, c1, c2, q,mpl::int_<PatternContext>()); \
            }                                                           \
            template<int PatternContext> \
                value_type                                              \
                evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<PatternContext>, mpl::int_<1>  ) const noexcept \
            {                                                           \
                value_type res( value_type( 0 ) );                      \
                for(uint16_type ii = 0; ii < l_type::shape::N; ++ii ) \
                    res += M_left.evalijq(i, j, c1, ii, q,mpl::int_<PatternContext>() ) VF_OP_SYMBOL( O ) M_right.evalijq(i, j, ii, c2, q,mpl::int_<PatternContext>() ); \
                return res;                                             \
            }                                                           \
                                               \
                value_type                                              \
                evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q  ) const noexcept \
            {                                                           \
                return evaliq(i, c1, c2, q, mpl::int_<shape_op>() );    \
            }                                                           \
                                               \
                value_type                                              \
                    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0>  ) const noexcept \
            {                                                           \
                if constexpr ( is_zero::value )                         \
                    return value_type( 0 );                             \
                else if constexpr ( !is_zero::update_and_eval_left && is_zero::update_and_eval_right ) \
                    return M_right.evaliq(i, c1, c2, q);               \
                else if constexpr ( is_zero::update_and_eval_left && !is_zero::update_and_eval_right ) \
                    return M_left.evaliq(i, c1, c2, q);                \
                else if constexpr ( l_type::shape::is_scalar && !r_type::shape::is_scalar ) \
                    return  M_left.evaliq(i, 0, 0, q) VF_OP_SYMBOL( O ) M_right.evaliq(i, c1, c2, q); \
                else if constexpr ( !l_type::shape::is_scalar && r_type::shape::is_scalar ) \
                    return  M_left.evaliq(i, c1, c2, q) VF_OP_SYMBOL( O ) M_right.evaliq(i, 0, 0, q); \
                else                                                    \
                    return  M_left.evaliq(i, c1, c2, q) VF_OP_SYMBOL( O ) M_right.evaliq(i, c1, c2, q); \
            }                                                           \
                                               \
                value_type                                              \
                evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1>  ) const noexcept \
            {                                                           \
                if constexpr ( is_zero::value )                         \
                    return value_type( 0 );                             \
                else                                                    \
                    {                                                   \
                        value_type res( value_type( 0 ) );              \
                        for(uint16_type ii = 0; ii < l_type::shape::N; ++ii ) \
                            res += M_left.evaliq(i, c1, ii, q ) VF_OP_SYMBOL( O ) M_right.evaliq(i, ii, c2, q ); \
                        return res;                                     \
                    }                                                   \
            }                                                           \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q ) const noexcept \
            {                                                           \
                return evalq( c1, c2, q, mpl::int_<shape_op>() );       \
            }                                                           \
                                                                        \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<0> ) const noexcept \
            {                                                           \
                if constexpr ( l_type::shape::is_scalar && !r_type::shape::is_scalar ) \
                    return M_left.evalq( 0, 0, q ) VF_OP_SYMBOL( O ) M_right.evalq( c1, c2, q ); \
                else if constexpr ( !l_type::shape::is_scalar && r_type::shape::is_scalar ) \
                    return M_left.evalq( c1, c2, q ) VF_OP_SYMBOL( O ) M_right.evalq( 0, 0, q ); \
                else                                                    \
                    return M_left.evalq( c1, c2, q ) VF_OP_SYMBOL( O ) M_right.evalq( c1, c2, q ); \
            }                                                           \
            value_type                                                  \
                evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<1> ) const noexcept \
            {                                                           \
                value_type res( value_type( 0 ) );                      \
                for(uint16_type ii = 0; ii < l_type::shape::N; ++ii )   \
                    res += M_left.evalq( c1, ii, q ) VF_OP_SYMBOL( O ) M_right.evalq( ii, c2, q ); \
                return res;                                             \
            }                                                           \
            l_type M_left;                                             \
            r_type M_right;                                            \
        }; /* tensor */                                                 \
        evaluate_type                                                   \
            evaluate(bool p,  worldcomm_ptr_t const& worldcomm ) const  \
        {                                                               \
            auto leval = M_left.evaluate(p,worldcomm).template cast<value_type>(); \
            auto reval = M_right.evaluate(p,worldcomm).template cast<value_type>(); \
            if ( leval.rows() == 1 && leval.cols() == 1 )               \
            {                                                           \
                if ( reval.rows() == 1 && reval.cols() == 1 )           \
                    return Eigen::Matrix<value_type,1,1>::Constant( leval(0,0) VF_OP_SYMBOL( O ) reval(0,0) ); \
                else if constexpr( is_op_mul )                          \
                    return leval(0,0) VF_OP_SYMBOL( O ) reval; \
                else                                                    \
                {                                                       \
                    CHECK( false ) << "should not go here";             \
                    return evaluate_type::Constant( 1,1,0. );           \
                }                                                       \
            }                                                           \
            else if ( reval.rows() == 1 && reval.cols() == 1 )          \
            {                                                           \
                if constexpr( is_op_mul || is_op_div )                  \
                    return leval VF_OP_SYMBOL( O ) reval(0,0); \
                else                                                    \
                {                                                       \
                    CHECK( false ) << "should not go here";             \
                    return evaluate_type::Constant( 1,1,0. );           \
                }                                                       \
            }                                                           \
            else if constexpr( (L_type::evaluate_type::SizeAtCompileTime == Eigen::Dynamic || \
                                R_type::evaluate_type::SizeAtCompileTime == Eigen::Dynamic) && !is_op_div) \
                                 return leval VF_OP_SYMBOL( O ) reval;  \
            else if constexpr( L_type::evaluate_type::SizeAtCompileTime > 1 && R_type::evaluate_type::SizeAtCompileTime > 1 ) \
            {                                                           \
                return leval VF_OP_SYMBOL( O ) reval;                   \
            }                                                           \
            else                                                        \
            {                                                           \
                CHECK( false ) << "should not go here";                 \
                return Eigen::Matrix<value_type,1,1>::Constant( 0 );    \
            }                                                           \
        }                                                               \
                                                                        \
        std::string expressionStr() const                               \
        {                                                               \
            return std::string();/*M_left.expressionStr() + BOOST_PP_STRINGIZE( VF_OP_SYMBOL( O ) ) + M_right.expressionStr();*/ \
        }                                                               \
        BOOST_PP_IF(1,                                                  \
                    VF_SYMETRIC,                                        \
                    BOOST_PP_EMPTY )()                                  \
            BOOST_PP_IF(BOOST_PP_AND(VF_OP_IS_ADD(O),                   \
                                     BOOST_PP_AND(VF_TYPE_IS_EXPR(L),VF_TYPE_IS_EXPR(R))), \
                        VF_ASSEMBLE,                                    \
                        BOOST_PP_EMPTY )()                              \
                                                                        \
            protected:                                                  \
            VF_OP_NAME( O )() {}                                        \
                                                                        \
        L_type M_left;                                                  \
        R_type M_right;                                                \
    };                                                                  \
    template <BOOST_PP_IF( VF_TYPE_IS_EXPR( L ),                        \
                           BOOST_PP_IDENTITY(class VF_TYPE_NAME(L)),    \
                           BOOST_PP_EMPTY                               \
                           )()                                          \
              BOOST_PP_IF( BOOST_PP_AND( VF_TYPE_IS_EXPR(L),            \
                                         VF_TYPE_IS_EXPR(R) ),          \
                           BOOST_PP_COMMA,                              \
                           BOOST_PP_EMPTY )()                           \
              BOOST_PP_IF( VF_TYPE_IS_EXPR( R ),                        \
                           BOOST_PP_IDENTITY(class VF_TYPE_NAME(R)),    \
                           BOOST_PP_EMPTY                               \
                           )()>                                         \
    inline                                                              \
        Expr< VF_OP_NAME( O )< VF_TYPE_NAME(L), VF_TYPE_NAME(R) > >     \
    operator VF_OP_SYMBOL( O )( VF_TYPE_TYPE_EXPR(L) VF_TYPE_CV(L) v, VF_TYPE_TYPE_EXPR(R) VF_TYPE_CV(R) w ) \
    {                                                                   \
        typedef VF_OP_NAME( O )<VF_TYPE_NAME( L ), VF_TYPE_NAME( R )> expr_t; \
        return Expr<expr_t> (expr_t ( VF_TYPE_TYPE_CST(L)(v) , VF_TYPE_TYPE_CST(R)(w) )); \
    }                                                                   \
    /**/
#
# define VF_ASSEMBLE()                                                  \
    template<typename Elem1, typename Elem2, typename FormType>         \
    void assemble( std::shared_ptr<Elem1> const& __u,  std::shared_ptr<Elem2> const& __v, FormType& __f ) const \
    {                                                                   \
        M_left.assemble( __u, __v, __f );                              \
        M_right.assemble( __u, __v, __f );                             \
    }                                                                   \
                                                                        \
    template<typename Elem1, typename FormType>                         \
    void assemble( std::shared_ptr<Elem1> const& __v, FormType& __f  ) const \
    {                                                                   \
        M_left.assemble( __v, __f );                                   \
        M_right.assemble( __v, __f );                                  \
    }                                                                   \
    /**/
#
# define VF_SYMETRIC()                                              \
    bool isSymetric() const                                         \
        {                                                           \
            return false;                                           \
        }                                                           \
 /**/
#

namespace Feel
{
namespace vf
{

BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_ARRAY_OP_DECLARATION, 3, ( VF_APPLICATIVE_BINARY_OPS, VF_EXPRL_TYPES, VF_EXPRR_TYPES ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_ARRAY_OP_DECLARATION, 3, ( VF_APPLICATIVE_BINARY_OPS, VF_EXPRL_TYPES, VF_BUILTIN_TYPES ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_ARRAY_OP_DECLARATION, 3, ( VF_APPLICATIVE_BINARY_OPS, VF_BUILTIN_TYPES, VF_EXPRR_TYPES ) )

BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_ARRAY_OP, 3, ( VF_APPLICATIVE_BINARY_OPS, VF_EXPRL_TYPES, VF_EXPRR_TYPES ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_ARRAY_OP, 3, ( VF_APPLICATIVE_BINARY_OPS, VF_EXPRL_TYPES, VF_BUILTIN_TYPES ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( VF_BINARY_ARRAY_OP, 3, ( VF_APPLICATIVE_BINARY_OPS, VF_BUILTIN_TYPES, VF_EXPRR_TYPES ) )
}
}

/// \endcond detail
#endif /* FEELPP_OPERATIONS_HPP */
