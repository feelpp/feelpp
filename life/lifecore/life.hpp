/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-12-30

  Copyright (C) 2006,2007,2008,2009 Université de Grenoble 1

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
   \file life.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-12-30
 */
# ifndef __cplusplus
# error You must use C++ for Life
# endif

# ifndef _LIFE_HH_
# define _LIFE_HH_

#include <complex>

#include <boost/mpl/multiplies.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/lower_bound.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/sizeof.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/base.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/comparison.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


#include <boost/program_options.hpp>

#include <boost/cstdint.hpp>


#include <cmath>
#include <numeric>
#include <string>
#include <limits>
#include <iosfwd>

#include <lifeconfig.h>
#include <life/lifecore/info.hpp>
#include <life/lifecore/lifemacros.hpp>
#include <life/lifecore/lifeassert.hpp>

#include <life/lifecore/flags.hpp>

namespace Life
{
namespace mpl = boost::mpl;
namespace lambda = boost::lambda;
namespace po = boost::program_options;


namespace detail
{
/*
  "inline" is used for ignore_unused_variable_warning()
   and function_requires() to make sure there is no
   overhead with g++.
   (ripped from boost/concept_check.hpp in order to avoid importing a huge header)

   \code
   void f( int q )
   {
     // q is not used before the check
     LIFE_ASSERT( q ).error( "q <= 0" );
     // q is not used after the check

     // this ensures that in -DNDEBUG there are no warnings about
     // unused parameter
     detail::ignore_unused_variable_warning( q );
   }
   \endcode
 */
template <class T> inline void ignore_unused_variable_warning(const T&) { }
};

/*!  \page types_page Life Types
  \section types Types
  \subsection real Real Numbers

  Life defines a number of types that are used in the library.

  -# \c Real 64 bits real number type

  \section ints Integers

  Life defines a number of integer type that have controlled bit
  size. These types are constructed automatically by Life in order to have
  platform independant integer types.

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

  Life defines a number of useful aliases for integers
  -# \c dim_type an alias to uint16_type used to identify dimensions
  -# \c size_type an alias to size_t used as indices for arrays, vectors or matrices

*/

typedef double Real;
typedef double scalar_type;
typedef std::complex<double> complex_type;

//
// Create type that are machine independent
// \warning should test here the boost version
//
// only available in boost 1.32
#include <boost/mpl/eval_if.hpp>

/*! \namespace detail
  \internal
*/
namespace detail
{
template<int bit_size>
class no_int
{
private:
    no_int();
};

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
}
#if 0
typedef detail::integer<1>::type  int1_type;

typedef detail::integer<8>::type  int8_type;
typedef detail::integer<16>::type int16_type;
typedef detail::integer<32>::type int32_type;
typedef detail::integer<64>::type int64_type;
typedef detail::integer<128>::type int128_type;
#else
typedef boost::int8_t   int8_type;
typedef boost::int16_t  int16_type;
typedef boost::int32_t  int32_type;
#if !defined( BOOST_NO_INT64_T )
    typedef boost::int64_t  int64_type;
#endif // BOOST_NO_INT64_T
#endif // 0
typedef detail::real<32>::type real32_type;
typedef detail::real<64>::type real64_type;
typedef detail::real<96>::type real96_type;

BOOST_STATIC_ASSERT( ( boost::is_same<real32_type, float>::value ) );
BOOST_STATIC_ASSERT( ( boost::is_same<real64_type, double>::value ) );

// don't check for this, it fails on armel arch.
#if 0
BOOST_STATIC_ASSERT( ( boost::is_same<real96_type, long double>::value ) );
#endif // 0

/*! \namespace detail
  \internal
*/
namespace detail
{
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
}
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
typedef size_t size_type;
//BOOST_MPL_ASSERT_MSG( ( boost::is_same<size_type,uint32_type>::value), NOT_SAME, (size_t, long unsigned int, uint32_type, uint64_type));
/**
 * @name Constants
 */
//@{

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
 * Invalid size type value
 */
const size_type invalid_size_type_value = size_type( -1 );

//@}

} // end namespace Life

#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Life
{
// alias for program_options namespace
namespace po = boost::program_options;

// alias for date_time namespaces
namespace posix_time = boost::posix_time;
namespace gregorian = boost::gregorian;
}

#if defined( HAVE_ARPREC)
# define LIFE_HAVE_MP_REAL 1
#include <mp/mpreal.h>
#include <mp/mp.h>
#endif /* HAVE_ARPREC */


#if defined( HAVE_QDLIB ) || defined( HAVE_QD_H )
# define LIFE_HAVE_DD_REAL 1
# define LIFE_HAVE_QD_REAL 1
# include <qd/dd.h>
# include <qd/qd.h>
# include <qd/fpu.h>


/// numeric_limits<dd_real> specialization.
namespace std {
template<>
struct numeric_limits<dd_real>
{
      static const bool is_specialized = true;

      static const int digits = 32;
      static const int digits10 = 32;
      static const bool is_signed = true;
      static const bool is_integer = false;
      static const bool is_exact = false;
      static const int radix = __FLT_RADIX__;
      static dd_real epsilon() throw()
      { return dd_real::_eps; }
      static dd_real round_error() throw()
      { return 0.5; }
      static dd_real min() throw()
      { return dd_real::_eps; }
    };

/// numeric_limits<dd_real> specialization.
template<>
struct numeric_limits<qd_real>
{
      static const bool is_specialized = true;

      static const int digits = 64;
      static const int digits10 = 64;
      static const bool is_signed = true;
      static const bool is_integer = false;
      static const bool is_exact = false;
      static const int radix = __FLT_RADIX__;
      static qd_real epsilon() throw()
      { return qd_real::_eps; }
      static qd_real round_error() throw()
      { return 0.5; }

      static qd_real min() throw()
      { return qd_real::_eps; }


    };
}

#endif /* HAVE_QD */

#if defined( HAVE_MPFR )
#include <life/lifecore/mpfr.hpp>

namespace Life
{
/**
 *\brief C++ wrapper for multiple precision floating type from mpfr
 */
typedef mpfr::MpfrClass mp_type;

/**
 *\brief precision type for multiple-precision floating type
 *
 *\
 */
typedef mp_prec_t mp_precision_type;

/**
 *\brief set the precision for multiple precision floating-point type
 * \param prec required precision of computation
 *
 * \code
 * mp_precision_type prec = 200;
 * setMpPrecision (prec);
 * \endcode
 */
inline void setMpPrecision( mp_precision_type __prec )
{
    mpfr_set_default_prec (__prec);
}

const mp_type mp_eps = mpfr::pow( mp_type(  2 ), -mp_type::GetDefaultPrecision()+1 );

}
#endif // HAVE_MPFR


#include <life/lifecore/debug.hpp>

#include <boost/shared_ptr.hpp>


#if defined(HAVE_OPENMP)
#include <omp.h>

#define OMP_SET_NUM_THREADS(num) omp_set_num_threads(num)
#define OMP_GET_NUM_THREADS      omp_get_num_threads()
#define OMP_SET_DYNAMIC(num)     omp_set_dynamic(num)
#define OMP_SET_NESTED(num)      omp_set_nested(num)

// openmp timing
#define OMP_GET_WTIME            omp_get_wtime()
#define OMP_GET_WTICK            omp_get_wtick()

#else
#define OMP_SET_NUM_THREADS(num)
#define OMP_GET_NUM_THREADS     1
#define OMP_SET_DYNAMIC(num)
#define OMP_SET_NESTED(num)

#define OMP_GET_WTIME           0
#define OMP_GET_WTICK           0

#endif /* HAVE_OPENMP */



#endif

