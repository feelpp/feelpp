/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-08-17

  Copyright (C) 2005,2006 EPFL

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
   \file life/lifealg/traits.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-08-17
 */
#ifndef __LIFEALG_TRAITS_HPP
#define __LIFEALG_TRAITS_HPP 1

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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>


#include <life/lifecore/life.hpp>
#include <life/lifecore/traits.hpp>


namespace Life
{
namespace glas
{
namespace ublas = boost::numeric::ublas;

struct vector_tag {};
struct matrix_tag {};

/**
 * base template traits that serves also as fallback if no
 * specialization has been found. This might not work properly of
 * course (ie may fail at compile time). To remedy this issue one must
 * add the new type to the list of specialized vector or matrix that
 * are supported by Life.
 */
template<typename T>
struct traits
{
    typedef T self_type;
    typedef typename self_type::value_type value_type;
};

/*
 * Vector
 */
# define LIFE_GLAS_TRAITS_VECTOR_TYPE(T)              BOOST_PP_TUPLE_ELEM(5, 0 , T)
# define LIFE_GLAS_TRAITS_VECTOR_ITERATOR(T)          BOOST_PP_TUPLE_ELEM(5, 1 , T)
# define LIFE_GLAS_TRAITS_VECTOR_CONST_ITERATOR(T)    BOOST_PP_TUPLE_ELEM(5, 2 , T)
# define LIFE_GLAS_TRAITS_VECTOR_SIZE(T)              BOOST_PP_TUPLE_ELEM(5, 3 , T)
# define LIFE_GLAS_TRAITS_VECTOR_RESIZE(T)            BOOST_PP_TUPLE_ELEM(5, 4 , T)
#
# define LIFE_GLAS_TRAITS_VECTOR_TYPES                                 \
  BOOST_PP_TUPLE_TO_LIST(                                               \
      2,                                                                \
      (                                                                 \
          ( ublas::vector, iterator, const_iterator, size, resize  ),   \
          ( std::vector  , iterator, const_iterator, size, resize )     \
      )                                                                 \
  )                                                                     \
  /**/
#
# /* Generates code for all integral types. */
# define LIFE_GLAS_TRAITS_VECTOR_OP(_, T) \
      LIFE_GLAS_TRAITS_VECTOR_OP_CODE T   \
   /**/
#
#define LIFE_GLAS_TRAITS_VECTOR_OP_CODE(T,V)                                                   \
template<>                                                                                      \
struct traits<LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> >                      \
{                                                                                               \
    typedef LIFE_TRAITS_TYPE( T ) value_type;                                                  \
    typedef LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> self_type;               \
    typedef self_type::LIFE_GLAS_TRAITS_VECTOR_ITERATOR( V ) iterator;                         \
    typedef  self_type::LIFE_GLAS_TRAITS_VECTOR_CONST_ITERATOR( V ) const_iterator;            \
    typedef vector_tag type_tag;                                                                \
    static const bool is_vector = true;                                                         \
    static const bool is_matrix = false;                                                        \
                                                                                                \
};                                                                                              \
traits<LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> >::const_iterator             \
begin( LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> const& t )                    \
{                                                                                               \
    return t.begin();                                                                           \
}                                                                                               \
traits<LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> >::iterator                   \
begin( LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )>& t )                          \
{                                                                                               \
    return t.begin();                                                                           \
}                                                                                               \
traits<LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> >::const_iterator             \
end( LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> const& t )                      \
{                                                                                               \
    return t.end();                                                                             \
}                                                                                               \
traits<LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )> >::iterator                   \
end( LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )>& t )                            \
{                                                                                               \
    return t.end();                                                                             \
}                                                                                               \
inline                                                                                          \
size_type size( LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )>  const& t )          \
{                                                                                               \
    return t.LIFE_GLAS_TRAITS_VECTOR_SIZE(V)();                                                \
}                                                                                               \
inline                                                                                          \
void resize( LIFE_GLAS_TRAITS_VECTOR_TYPE( V )<LIFE_TRAITS_TYPE( T )>& t, size_type newsize ) \
{                                                                                               \
    return t.LIFE_GLAS_TRAITS_VECTOR_RESIZE(V)( newsize );                                     \
}                                                                                               \
 /**/


/**
 * Generate the type traits
 */
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_GLAS_TRAITS_VECTOR_OP, 2, (LIFE_TRAITS_TYPES, LIFE_GLAS_TRAITS_VECTOR_TYPES) );

/*
 * Matrices
 */
#
# define LIFE_GLAS_TRAITS_MATRIX_TYPE(T)              BOOST_PP_TUPLE_ELEM(5, 0 , T)
# define LIFE_GLAS_TRAITS_MATRIX_SIZE1(T)             BOOST_PP_TUPLE_ELEM(5, 1 , T)
# define LIFE_GLAS_TRAITS_MATRIX_SIZE2(T)             BOOST_PP_TUPLE_ELEM(5, 2 , T)
# define LIFE_GLAS_TRAITS_MATRIX_RESIZE(T)            BOOST_PP_TUPLE_ELEM(5, 3 , T)
# define LIFE_GLAS_TRAITS_MATRIX_NNZ(T)               BOOST_PP_TUPLE_ELEM(5, 4 , T)
#
#define UBLAS_MATRIX_ROW( T ) ublas::matrix<T, ublas::row_major>
#define UBLAS_MATRIX_COL( T ) ublas::matrix<T, ublas::column_major>
#define UBLAS_MATRIX_SPARSE_ROW( T ) ublas::compressed_matrix<T, ublas::row_major>
#define UBLAS_MATRIX_SPARSE_COL( T ) ublas::compressed_matrix<T, ublas::column_major>
#
# define LIFE_GLAS_TRAITS_MATRIX_TYPES                                         \
  BOOST_PP_TUPLE_TO_LIST(                                                       \
      4,                                                                        \
      (                                                                         \
          ( UBLAS_MATRIX_ROW       , size1, size2, resize, boost::none_t  ),    \
          ( UBLAS_MATRIX_COL       , size1, size2, resize, boost::none_t  ),    \
          ( UBLAS_MATRIX_SPARSE_ROW, size1, size2, resize, nnz  ),              \
          ( UBLAS_MATRIX_SPARSE_COL, size1, size2, resize, nnz  )               \
      )                                                                         \
  )                                                                             \
  /**/
#
# /* Generates code for all integral types. */
# define LIFE_GLAS_TRAITS_MATRIX_OP(_, T) \
      LIFE_GLAS_TRAITS_MATRIX_OP_CODE T   \
   /**/
#
#define LIFE_GLAS_TRAITS_MATRIX_OP_CODE(T,V)                                                                           \
template<>                                                                                                              \
struct traits<LIFE_GLAS_TRAITS_MATRIX_TYPE( V )( LIFE_TRAITS_TYPE( T ) ) >                                            \
{                                                                                                                       \
    typedef LIFE_TRAITS_TYPE( T ) value_type;                                                                          \
    typedef LIFE_GLAS_TRAITS_MATRIX_TYPE( V )( LIFE_TRAITS_TYPE( T ) ) self_type;                                     \
    typedef matrix_tag type_tag;                                                                                        \
    static const bool is_vector = false;                                                                                \
    static const bool is_matrix = true;                                                                                 \
};                                                                                                                      \
inline                                                                                                                  \
size_type nrows( LIFE_GLAS_TRAITS_MATRIX_TYPE( V )( LIFE_TRAITS_TYPE( T ) )  const& t )                               \
{                                                                                                                       \
    return t.LIFE_GLAS_TRAITS_MATRIX_SIZE1(V)();                                                                       \
}                                                                                                                       \
inline                                                                                                                  \
size_type ncols( LIFE_GLAS_TRAITS_MATRIX_TYPE( V )( LIFE_TRAITS_TYPE( T ) )  const& t )                               \
{                                                                                                                       \
    return t.LIFE_GLAS_TRAITS_MATRIX_SIZE2(V)();                                                                       \
}                                                                                                                       \
inline                                                                                                                  \
void resize( LIFE_GLAS_TRAITS_MATRIX_TYPE( V )( LIFE_TRAITS_TYPE( T ) )& t, size_type newsize1, size_type newsize2 )  \
{                                                                                                                       \
    return t.LIFE_GLAS_TRAITS_MATRIX_RESIZE(V)( newsize1, newsize2 );                                                  \
}                                                                                                                       \
 /**/


/**
 * Generate the type traits
 */
BOOST_PP_LIST_FOR_EACH_PRODUCT( LIFE_GLAS_TRAITS_MATRIX_OP, 2, (LIFE_TRAITS_TYPES, LIFE_GLAS_TRAITS_MATRIX_TYPES) );

} // glas

} // Life

#endif /*__LIFEALG_TRAITS_HPP */
