//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 17 Jun 2018
//! @copyright 2018-2022 Feel++ Consortium
//!

#ifndef FEELPP_FEELTYPES_HPP
#define FEELPP_FEELTYPES_HPP 1

#include <limits>
#include <boost/cstdint.hpp>


/**
 *\page Types Types
   \ingroup Core
   \section types Types
   \subsection real Real Numbers

  Feel defines a number of types that are used in the library.

  -# \c real32_type 32 bits real number type
  -# \c real64_type 64 bits real number type

  \section ints Integers

  Feel defines a number of integer type that have controlled bit
  size. These types are constructed automatically by Feel in order to have
  platform independent integer types.

  Here is the list of signed integers:

  -# \c int1_type  a 1 bit signed integer
  -# \c int8_type a 8 bit signed integer
  -# \c int16_type a 16 bit signed integer
  -# \c int32_type a 32 bit signed integer
  -# \c int64_type a 64 bit signed integer

  Here is the list of unsigned integers:

  -# \c uint1_type a 1 bit unsigned integer
  -# \c uint8_type a 8 bit unsigned integer
  -# \c uint16_type a 16 bit unsigned integer
  -# \c uint32_type a 32 bit unsigned integer
  -# \c uint64_type a 64 bit unsigned integer

  Feel defines a number of useful aliases for integers
  -# \c dim_type an alias to uint16_type used to identify dimensions
  -# \c rank_type an alias to uint16_type used for mpi process ranks
  -# \c size_type an alias to size_t used as indices for arrays, vectors or matrices

  Feel++ defines also some invalid number associated to each unsigned integer
  types. it usually serves to identify a situation where something went
  wrong. The invalid value corresponds to the maximum value of the unsigned type.
  -# invalid_uint8_type_value
  -# invalid_uint16_type_value
  -# invalid_uint32_type_value
  -# invalid_uint64_type_value
  -# invalid_dim_type_value
  -# invalid_rank_type_value
  -# invalid_v<size_type>

*/

typedef double Real;
typedef double float64_t;
typedef double scalar_type;

//
// Create type that are machine independent
// \warning should test here the boost version
//
// only available in boost 1.32
// #include <boost/mpl/eval_if.hpp>

namespace Feel
{
//! @cond
namespace detail
{
template<int bit_size>
class no_int
{
private:
    no_int();
};

#if 0
template< int bit_size >
struct integer
{
    typedef mpl::list<signed char,signed short, signed int, signed long int, signed long long> builtins_;
    typedef typename mpl::base< typename mpl::lower_bound<
    mpl::transform_view< builtins_, mpl::multiplies< mpl::sizeof_<mpl::placeholders::_1>, mpl::int_<8> >
    >
    , mpl::integral_c<size_t, bit_size>
    >::type >::type iter_;

    typedef typename mpl::end<builtins_>::type last_;
    typedef typename mpl::eval_if<
    boost::is_same<iter_,last_>
    , mpl::identity< no_int<bit_size> >
    , mpl::deref<iter_>
    >::type type;
};

template< int bit_size >
struct real
{
    typedef mpl::list<float, double, long double> builtins_;
    typedef typename mpl::base< typename mpl::lower_bound<
    mpl::transform_view< builtins_, mpl::multiplies< mpl::sizeof_<mpl::placeholders::_1>, mpl::int_<8> >
    >
    , mpl::integral_c<size_t, bit_size>
    >::type >::type iter_;

    typedef typename mpl::end<builtins_>::type last_;
    typedef typename mpl::eval_if<
    boost::is_same<iter_,last_>
    , mpl::identity< no_int<bit_size> >
    , mpl::deref<iter_>
    >::type type;
};
#endif
}
//! @endcond
#if 0
typedef detail::integer<1>::type  int1_type;

typedef detail::integer<8>::type  int8_type;
typedef detail::integer<16>::type int16_type;
typedef detail::integer<32>::type int32_type;
typedef detail::integer<64>::type int64_type;
typedef detail::integer<128>::type int128_type;

typedef detail::real<32>::type real32_type;
typedef detail::real<64>::type real64_type;

#else
typedef boost::int8_t   int8_type;
typedef boost::int16_t  int16_type;
typedef boost::int32_t  int32_type;

typedef float real32_type;
typedef double real64_type;

#if !defined( BOOST_NO_INT64_T )
typedef boost::int64_t  int64_type;
#endif // BOOST_NO_INT64_T
#endif // 0


//! @cond
namespace detail
{
#if 0
template< int bit_size >
struct unsigned_integer
{
    //typedef mpl::list<unsigned char,unsigned short, long unsigned int, long unsigned int,  long unsigned long> builtins_;
    typedef mpl::list<unsigned char,unsigned short, unsigned int, unsigned long int,  unsigned long long> builtins_;
    typedef typename mpl::base< typename mpl::lower_bound<
    mpl::transform_view< builtins_
    , mpl::multiplies< mpl::sizeof_<mpl::placeholders::_1>, mpl::int_<8> >
    >
    , mpl::integral_c<size_t, bit_size>
    >::type >::type iter_;

    typedef typename mpl::end<builtins_>::type last_;
    typedef typename mpl::eval_if<
    boost::is_same<iter_,last_>
    , mpl::identity< no_int<bit_size> >
    , mpl::deref<iter_>
    >::type type;
};
#endif
}
//! @endcond

#if 0
typedef detail::unsigned_integer<1>::type  uint1_type;
typedef detail::unsigned_integer<8>::type  uint8_type;
typedef detail::unsigned_integer<16>::type uint16_type;
typedef detail::unsigned_integer<32>::type uint32_type;
typedef detail::unsigned_integer<64>::type uint64_type;
typedef detail::unsigned_integer<128>::type uint128_type;
#else
typedef boost::uint8_t   uint8_type;
typedef boost::uint16_t  uint16_type;
typedef boost::uint32_t  uint32_type;
#if !defined( BOOST_NO_INT64_T )
typedef boost::uint64_t  uint64_type;

/**
 * @typedef int64_type marker_type
 * marker_type is the type used to store the geometric entity flags
 * It is a 64 bits integer so thet we have a lot of room
 */
typedef int64_type flag_type;

const int64_type invalid_flag_type_value = std::numeric_limits<int32_type>::min();
#else
typedef int32_type flag_type;

const int32_type invalid_flag_type_value = std::numeric_limits<int32_type>::min();
#endif // BOOST_NO_INT64_T
#endif
//! dimension type
typedef uint16_type dim_type;

//! Indices (starting from 0)
//typedef size_t size_type;
using size_type = size_t;
//using size_type = uint32_type;

using index_type = uint32_type;

//! dof id type
using dof_id_type = index_type;

//! type for mpi rank ids
typedef uint16_type rank_type;

//! quadrature order type
typedef uint16_type quad_order_type;

#if defined( __APPLE__ )
typedef unsigned int uint;
#endif // __APPLE__ç
//BOOST_MPL_ASSERT_MSG( ( boost::is_same<size_type,uint32_type>::value), NOT_SAME, (size_t, long unsigned int, uint32_type, uint64_type));
/**
 * @name Constants
 */
//@{

//!
//! get invalid value from type \p T
//!
template<typename T>
inline const T invalid_v = T(-1);

/**
 * Invalid uint8_type value
 */
const uint8_type invalid_uint8_type_value = uint8_type( -1 );

/**
 * Invalid uint16_type value
 */
const uint16_type invalid_uint16_type_value = uint16_type( -1 );

/**
 * Invalid uint32_type value
 */
const uint32_type invalid_uint32_type_value = uint32_type( -1 );

#if !defined( BOOST_NO_INT64_T )
/**
 * Invalid uint64_type value
 */
const uint64_type invalid_uint64_type_value = uint64_type( -1 );
#endif  // BOOST_NO_INT64_T

/**
 * Invalid dim type value
 */
const dim_type invalid_dim_type_value = dim_type( -1 );

/**
 * Invalid dim type value
 */
const rank_type invalid_rank_type_value = rank_type( -1 );


/**
 * Invalid size type value
 */
const size_type invalid_size_type_value = size_type( -1 );

/**
 * Quadrature order deduced from the expression to integrate
 */
const uint16_type quad_order_from_expression = invalid_uint16_type_value;

//@}

//!
//! get the underlying value of an enum class entry
//!
template <typename E>
constexpr auto to_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}
template< typename E , typename T>
constexpr inline typename std::enable_if< std::is_enum<E>::value &&
                                          std::is_integral<T>::value, E
                                          >::type 
to_enum( T value ) noexcept 
{
    return static_cast<E>( value );
}
} // end namespace Feel

#endif
