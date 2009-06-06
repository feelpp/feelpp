/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-07-28

  Copyright (C) 2005,2006 EPFL

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
   \file typetraits.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-07-28
 */
#ifndef __typetraits_HPP
#define __typetraits_HPP 1

#include <boost/concept_check.hpp>

#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <boost/preprocessor/logical/and.hpp>
#include <boost/preprocessor/control/if.hpp>
#include <boost/preprocessor/list/at.hpp>
#include <boost/preprocessor/list/cat.hpp>
#include <boost/preprocessor/list/for_each_product.hpp>
#include <boost/preprocessor/logical/or.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/facilities/identity.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/list/filter.hpp>

#if !defined( LIFE_TRAITS_HPP)
//#error life/lifecore/typetraits.hpp must not be used directly, use life/lifecore/traits.hpp instead
#endif

#if defined( HAVE_QD_H )
qd_real floor(const qd_real &a);
qd_real ceil(const qd_real &a);
qd_real npwr(const qd_real &a, int n);
#endif

#if defined HAVE_ARPREC
inline mp_real_temp floor( mp_real const& m ) {
    if ( m >0 ) return ::aint( m );
    return ::aint( m-1 );
}
inline mp_real_temp ceil( mp_real const& m ) {
    if ( m >0 )return ::aint( m+1 );
    return ::aint( m );
}
#endif // HAVE_ARPREC

namespace Life
{

template<typename T>
struct type_traits
{
};

# define LIFE_FUNC_NAME(T)          BOOST_PP_TUPLE_ELEM(3, 0 , T)
# define LIFE_FUNC_CALL(T)          BOOST_PP_TUPLE_ELEM(3, 1 , T)
# define LIFE_FUNC_NONS(T)          BOOST_PP_TUPLE_ELEM(3, 2 , T)

#define LIFE_STD_FUNCS                                        \
BOOST_PP_TUPLE_TO_LIST(                                        \
    19,                                                        \
    (                                                          \
        (abs       , std::abs  , abs  ),   /* absolute value*/        \
        (sqrt      , std::sqrt , sqrt ),  /* square root */           \
        (norm1     , std::abs  , abs  ),   /* norm1    */             \
        (norm2     , std::abs  , abs  ),   /* norm2    */             \
        (norm_inf  , std::abs  , abs  ),   /* norm_inf */             \
        (cos       , std::cos  , cos  ),                              \
        (sin       , std::sin  , sin  ),                              \
        (tan       , std::tan  , tan  ),                              \
        (acos      , std::acos , acos ),                              \
        (asin      , std::asin , asin ),                              \
        (atan      , std::atan , atan ),                              \
        (cosh      , std::cosh , cosh ),                              \
        (sinh      , std::sinh , sinh ),                              \
        (tanh      , std::tanh , tanh ),                              \
        (exp       , std::exp  , exp  ),                              \
        (log       , std::log  , log  ),                              \
        (log10     , std::log10, log10),                              \
        (ceil      , std::ceil , ceil ),                              \
        (floor     , std::floor, floor)                               \
     )                                                         \
    )                                                          \
    /**/
#define LIFE_STDCOMPLEX_FUNCS                                  \
    BOOST_PP_TUPLE_TO_LIST(                                    \
                           11,                                        \
                           (                                            \
                            (abs       , std::abs  , abs  ),   /* absolute value*/ \
                            (sqrt      , std::sqrt , sqrt ),  /* square root */ \
                            (norm1     , std::abs  , abs  ),   /* norm1    */ \
                            (norm2     , std::abs  , abs  ),   /* norm2    */ \
                            (norm_inf  , std::abs  , abs  ),   /* norm_inf */ \
                            (cos       , std::cos  , cos  ),            \
                            (sin       , std::sin  , sin  ),            \
                            (tan       , std::tan  , tan  ),            \
                            (exp       , std::exp  , exp  ),            \
                            (log       , std::log  , log  ),            \
                            (log10     , std::log10, log10 )            \
                                                                        ) \
                                                               )        \
    /**/
#define LIFE_STD_BINARY_FUNCS                                  \
BOOST_PP_TUPLE_TO_LIST(                                        \
    1,                                                         \
    (                                                          \
        (pow       , std::pow, pow  )                          \
     )                                                         \
    )                                                          \
    /**/
#
#define LIFE_MP_FUNCS                                     \
BOOST_PP_TUPLE_TO_LIST(                                        \
    18,                                                        \
    (                                                          \
        (abs       , mpfr::abs  , abs  ),   /* absolute value*/           \
        (sqrt      , mpfr::sqrt , sqrt ),  /* square root */              \
        (norm1     , mpfr::abs  , abs  ),   /* norm1    */                \
        (norm2     , mpfr::abs  , abs  ),   /* norm2    */                \
        (norm_inf  , mpfr::abs  , abs  ),   /* norm_inf */                \
        (cos       , mpfr::cos  , cos  ),                                 \
        (sin       , mpfr::sin  , sin  ),                                 \
        (tan       , mpfr::tan  , tan  ),                                 \
        (acos      , mpfr::acos , acos ),                                 \
        (asin      , mpfr::asin , asin ),                                 \
        (atan      , mpfr::atan , atan ),                                 \
        (cosh      , mpfr::cosh , cosh ),                                 \
        (sinh      , mpfr::sinh , sinh ),                                 \
        (tanh      , mpfr::tanh , tanh ),                                 \
        (exp       , mpfr::exp  , exp  ),                                 \
        (log       , mpfr::log  , log  ),                                 \
        (ceil      , mpfr::ceil , ceil ),                                 \
        (floor     , mpfr::floor, floor)                                  \
     )                                                         \
    )                                                          \
    /**/
#
#
#define LIFE_MP_BINARY_FUNCS                                  \
BOOST_PP_TUPLE_TO_LIST(                                        \
    1,                                                         \
    (                                                          \
        (pow       , mpfr::pow, pow  )                           \
     )                                                         \
    )                                                          \
    /**/
#
#define LIFE_GLOBAL_FUNCS                                     \
BOOST_PP_TUPLE_TO_LIST(                                        \
    18,                                                        \
    (                                                          \
        (abs       , ::abs  , abs  ),   /* absolute value*/           \
        (sqrt      , ::sqrt , sqrt ),  /* square root */              \
        (norm1     , ::abs  , abs  ),   /* norm1    */                \
        (norm2     , ::abs  , abs  ),   /* norm2    */                \
        (norm_inf  , ::abs  , abs  ),   /* norm_inf */                \
        (cos       , ::cos  , cos  ),                                 \
        (sin       , ::sin  , sin  ),                                 \
        (tan       , ::tan  , tan  ),                                 \
        (acos      , ::acos , acos ),                                 \
        (asin      , ::asin , asin ),                                 \
        (atan      , ::atan , atan ),                                 \
        (cosh      , ::cosh , cosh ),                                 \
        (sinh      , ::sinh , sinh ),                                 \
        (tanh      , ::tanh , tanh ),                                 \
        (exp       , ::exp  , exp  ),                                 \
        (log       , ::log  , log  ),                                 \
        (ceil      , ::ceil , ceil ),                                 \
        (floor     , ::floor, floor)                                  \
     )                                                         \
    )                                                          \
    /**/
#
#define LIFE_GLOBAL_BINARY_FUNCS                              \
BOOST_PP_TUPLE_TO_LIST(                                        \
    1,                                                         \
    (                                                          \
        (pow       , ::npwr, npwr  )                           \
     )                                                         \
    )                                                          \
    /**/

#
#if 1
# define LIFE_TRAITS_FUNC_REAL(T, t)           \
  BOOST_PP_IF( LIFE_TRAITS_IS_COMPLEX(T),     \
               t.real(),  \
               t )      \
 /**/
#
# define LIFE_TRAITS_FUNC_IMAG(T, t)                                       \
  BOOST_PP_IF( LIFE_TRAITS_IS_COMPLEX(T),                                 \
               BOOST_PP_IDENTITY( t.imag() ),                              \
               BOOST_PP_IDENTITY( LIFE_TRAITS_REAL_TYPE(T)( 0.0 ) ) )()   \
 /**/
#
# define LIFE_TRAITS_FUNC_CONJ(T, t)                   \
  BOOST_PP_IF( LIFE_TRAITS_IS_COMPLEX(T),             \
               BOOST_PP_IDENTITY( std::conj( t ) ),    \
               BOOST_PP_IDENTITY( t ) )()              \
 /**/
#else
# define LIFE_TRAITS_FUNC_REAL(T, t) t
# define LIFE_TRAITS_FUNC_IMAG(T, t) 0.0
# define LIFE_TRAITS_FUNC_CONJ(T, t) t
#endif
#
# define LIFE_TRAITS_TYPE(T)              BOOST_PP_TUPLE_ELEM(7, 0 , T)
# define LIFE_TRAITS_REAL_TYPE(T)         BOOST_PP_TUPLE_ELEM(7, 1 , T)
# define LIFE_TRAITS_IS_FLOATING(T)       BOOST_PP_TUPLE_ELEM(7, 2 , T)
# define LIFE_TRAITS_IS_COMPLEX(T)        BOOST_PP_TUPLE_ELEM(7, 3 , T)
# define LIFE_TRAITS_RANK(T)              BOOST_PP_TUPLE_ELEM(7, 4 , T)
# define LIFE_TRAITS_EPSILON(T)           BOOST_PP_TUPLE_ELEM(7, 5 , T)
# define LIFE_TRAITS_FUNC_TYPE(T)         BOOST_PP_TUPLE_ELEM(7, 6 , T)
#
const double factor_from_eps = 50;
const float  factor_from_eps_fl = 50;

#if defined( HAVE_QD_H )
#define QD_DD_TYPE  ( dd_real, dd_real, 1, 0, 11, dd_real::_eps*factor_from_eps , 2 ),
#define QD_QD_TYPE  ( qd_real, qd_real, 1, 0, 12, qd_real::_eps*factor_from_eps , 2 ),
#define QD_NTYPES 2
#else
#define QD_NTYPES 0
#define QD_DD_TYPE
#define QD_QD_TYPE
#endif // HAVE_QD_H

#if defined( HAVE_MPFR )
#define MPFR_MP_TYPE ( mp_type, mp_type, 1, 0, 10 , mp_eps*factor_from_eps, 1 ),
#define MPFR_NTYPES 1
#else
#define MPFR_NTYPES 0
#define MPFR_MP_TYPE
#endif // HAVE_MPFR

#define LIFE_NUMERICAL_NTYPES BOOST_PP_ADD(6, BOOST_PP_ADD( QD_NTYPES, MPFR_NTYPES ) )
# define LIFE_TRAITS_TYPES \
   BOOST_PP_TUPLE_TO_LIST( \
      LIFE_NUMERICAL_NTYPES, \
      ( \
       QD_DD_TYPE                               \
       QD_QD_TYPE                               \
       MPFR_MP_TYPE                                                     \
       ( float      , float      , 1, 0, 7 , std::numeric_limits<float>::epsilon()*factor_from_eps_fl , 0 ), \
       ( double     , double     , 1, 0, 8 , std::numeric_limits<double>::epsilon()*factor_from_eps, 0 ), \
       ( real96_type, real96_type, 1, 0, 9 , std::numeric_limits<real96_type>::epsilon()*factor_from_eps, 0 ), \
       ( std::complex<float>, float, 1, 1, 10, std::numeric_limits<float>::epsilon()*factor_from_eps_fl , 3 ), \
       ( std::complex<double>, double, 1, 1, 11, std::numeric_limits<double>::epsilon()*factor_from_eps , 3 ), \
       ( std::complex<real96_type>, real96_type, 1, 1, 12, std::numeric_limits<real96_type>::epsilon()*factor_from_eps , 3 ) \
      ) \
   ) \
   /**/
#
# /* Generates code for all integral types. */
# define LIFE_TRAITS_OP(_, T) \
      LIFE_TRAITS_OP_CODE T   \
   /**/
#
#define LIFE_TRAITS_OP_CODE(T)                                                     \
template<>                                                                          \
struct type_traits<LIFE_TRAITS_TYPE( T )> {                                        \
    typedef type_traits<LIFE_TRAITS_TYPE( T )> self_type;                          \
    typedef LIFE_TRAITS_TYPE( T ) value_type;                                      \
    typedef const value_type &const_reference;                                      \
    typedef value_type &reference;                                                  \
    typedef LIFE_TRAITS_REAL_TYPE( T ) real_type;                                  \
    typedef LIFE_TRAITS_REAL_TYPE( T ) precision_type;                             \
                                                                                    \
    static const bool is_floating = LIFE_TRAITS_IS_FLOATING( T );       \
    static const uint16_type rank = LIFE_TRAITS_RANK( T );              \
    static real_type epsilon() { return LIFE_TRAITS_EPSILON( T ); }     \
};                                                                      \
inline                                                                              \
LIFE_TRAITS_REAL_TYPE( T ) real( LIFE_TRAITS_TYPE( T ) const& t )                 \
{                                                                                   \
    boost::ignore_unused_variable_warning(t);                                       \
    return LIFE_TRAITS_FUNC_REAL( T, t );                                          \
}                                                                                   \
inline                                                                              \
LIFE_TRAITS_REAL_TYPE( T ) imag ( LIFE_TRAITS_TYPE( T ) const& t )                \
{                                                                                   \
    boost::ignore_unused_variable_warning(t);                                       \
    return LIFE_TRAITS_FUNC_IMAG( T, t );                                          \
}                                                                                   \
inline                                                                              \
LIFE_TRAITS_TYPE( T ) conj ( LIFE_TRAITS_TYPE( T ) const& t)                      \
{                                                                                   \
    boost::ignore_unused_variable_warning(t);                                       \
    return LIFE_TRAITS_FUNC_CONJ( T, t );                                          \
}                                                                                   \
 /**/


/**
 * Generate the type traits
 */
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_TRAITS_OP, 1, (LIFE_TRAITS_TYPES) );

namespace math
{
#
# define LIFE_PRED_FUNC(d, data, elem) BOOST_PP_EQUAL(LIFE_TRAITS_FUNC_TYPE(elem), data)
#
/**
 * Generate the unary functions for each type
 */
# define LIFE_UNARY_FUNCS_OP(_,TF) \
  LIFE_UNARY_FUNCS_OP_CODE TF      \
  /**/
#
# define LIFE_UNARY_FUNCS_OP_CODE(T,F)                                             \
  inline                                                                            \
  LIFE_TRAITS_TYPE( T ) LIFE_FUNC_NAME( F )( LIFE_TRAITS_TYPE( T ) const& t )    \
  {                                                                                 \
      return LIFE_FUNC_CALL( F )( t );                                             \
  }                                                                                 \
  /**/
#
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_UNARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 0, LIFE_TRAITS_TYPES ), LIFE_STD_FUNCS) );
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_UNARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 1, LIFE_TRAITS_TYPES ), LIFE_MP_FUNCS) );
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_UNARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 2, LIFE_TRAITS_TYPES ), LIFE_GLOBAL_FUNCS) );
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_UNARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 3, LIFE_TRAITS_TYPES ), LIFE_STDCOMPLEX_FUNCS) );

/**
 * Generate the binary functions for each type
 */
# define LIFE_BINARY_FUNCS_OP(_,TF) \
  LIFE_BINARY_FUNCS_OP_CODE TF      \
  /**/
#
# define LIFE_BINARY_FUNCS_OP_CODE(T,F)                                                         \
  template<typename Type2>                                                                       \
  inline                                                                                         \
  LIFE_TRAITS_TYPE( T ) LIFE_FUNC_NAME( F )( LIFE_TRAITS_TYPE( T ) const& x, Type2 const& y ) \
  {                                                                                              \
      return LIFE_FUNC_CALL( F )( x, y );                                                       \
  }                                                                                              \
  /**/
#
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_BINARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 0, LIFE_TRAITS_TYPES ), LIFE_STD_BINARY_FUNCS) );
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_BINARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 1, LIFE_TRAITS_TYPES ), LIFE_MP_BINARY_FUNCS) );
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_BINARY_FUNCS_OP, 2, (BOOST_PP_LIST_FILTER(LIFE_PRED_FUNC, 2, LIFE_TRAITS_TYPES ), LIFE_GLOBAL_BINARY_FUNCS) );
} // math



} // Life

#if defined( HAVE_QD_H )
//
// make dd_real/qd_real known to boost::lambda
//
namespace boost {
namespace lambda {
namespace detail {

template <> struct promote_code<dd_real> { static const int value = 6000; };
template <> struct promote_code<qd_real> { static const int value = 7000; };

} // namespace detail
} // namespace lambda
} // namespace boost

#endif /* HAVE_QD_H */

#endif /* __typetraits_H */
