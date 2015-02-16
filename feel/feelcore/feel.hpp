/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2006-12-30

 Copyright (C) 2006,2007,2008,2009,2010 Universite de Grenoble 1
 Copyright (C) 2011-2015 Feel++ Consortium

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
 \file feel.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-12-30
 */
# ifndef __cplusplus
# error You must use C++ for Feel
# endif

# ifndef _FEELPP_HH_
# define _FEELPP_HH_



#if defined(__APPLE__)
#undef tolower
#undef toupper
#endif

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

#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <boost/assign/list_of.hpp>

#include <boost/math/constants/constants.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdivision-by-zero"
#endif
#include <boost/mpi.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <boost/program_options.hpp>

#include <boost/cstdint.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include <cmath>
#include <numeric>
#include <string>
#include <limits>
#include <iosfwd>

#if defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning(disable:780)
#endif
#include <glog/logging.h>
#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif

#include <feel/feelconfig.h>
#include <feel/feelcore/info.hpp>
#include <feel/feelcore/feelmacros.hpp>
#include <feel/feelcore/feelassert.hpp>

#include <feel/feelcore/flags.hpp>
#include <feel/feelcore/serialization.hpp>

#if defined( FEELPP_HAS_TBB )
#include <tbb/tick_count.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/mutex.h>
#endif // FEELPP_HAS_TBB



namespace Feel
{
namespace assign = boost::assign;
namespace fs = boost::filesystem;
namespace mpl = boost::mpl;
namespace lambda = boost::lambda;
namespace po = boost::program_options;

// bring boost.mpi into Feel realm
namespace mpi=boost::mpi;

namespace constants=boost::math::constants;
//namespace constants = boost::math::constants;
//using boost::math::double_constants;
const double pi = constants::pi<double>();
const double two_pi = constants::two_pi<double>();

namespace algorithm=boost::algorithm;
using google::WARNING;
using google::ERROR;
using google::INFO;
using google::FATAL;

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
     FEELPP_ASSERT( q ).error( "q <= 0" );
     // q is not used after the check

     // this ensures that in -DNDEBUG there are no warnings about
     // unused parameter
     detail::ignore_unused_variable_warning( q );
   }
   \endcode
 */
template <class T> inline void ignore_unused_variable_warning( const T& ) { }
}

/*!  \page Types Feel Types
  \section types Types
  \subsection real Real Numbers

  Feel defines a number of types that are used in the library.

  -# \c real32_type 32 bits real number type
  -# \c real64_type 64 bits real number type

  \section ints Integers

  Feel defines a number of integer type that have controlled bit
  size. These types are constructed automatically by Feel in order to have
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
  -# invalid_size_type_value

*/

typedef double Real;
typedef double float64_t;
typedef double scalar_type;

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

BOOST_STATIC_ASSERT( ( boost::is_same<real32_type, float>::value ) );
BOOST_STATIC_ASSERT( ( boost::is_same<real64_type, double>::value ) );

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

//! type for mpi rank ids
typedef uint16_type rank_type;

#if defined( __APPLE__ )
typedef unsigned int uint;
#endif // __APPLE__
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
 * Invalid dim type value
 */
const rank_type invalid_rank_type_value = rank_type( -1 );


/**
 * Invalid size type value
 */
const size_type invalid_size_type_value = size_type( -1 );

//@}

} // end namespace Feel

#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Feel
{
// alias for program_options namespace
namespace po = boost::program_options;

std::string
prefixvm( std::string const& prefix,
          std::string const& opt,
          std::string const& sep="." );


// alias for date_time namespaces
namespace posix_time = boost::posix_time;
namespace gregorian = boost::gregorian;

namespace meta
{
template<typename TheArgs>
struct remove_all
{
    typedef typename boost::remove_pointer<
    typename boost::remove_const<
    typename boost::remove_reference<
    TheArgs
    >::type
    >::type
    >::type type;
};
}
}

#if defined( FEELPP_HAS_ARPREC)
# define FEELPP_HAS_MP_REAL 1
#include <mp/mpreal.h>
#include <mp/mp.h>
#endif /* FEELPP_HAS_ARPREC */


#if defined( FEELPP_HAS_QDLIB ) || defined( FEELPP_HAS_QD_H )
# define FEELPP_HAS_DD_REAL 1
# define FEELPP_HAS_QD_REAL 1
# include <qd/dd.h>
# include <qd/qd.h>
# include <qd/fpu.h>


/// numeric_limits<dd_real> specialization.
namespace std
{
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
    {
        return dd_real::_eps;
    }
    static dd_real round_error() throw()
    {
        return 0.5;
    }
    static dd_real min() throw()
    {
        return dd_real::_eps;
    }
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
    {
        return qd_real::_eps;
    }
    static qd_real round_error() throw()
    {
        return 0.5;
    }

    static qd_real min() throw()
    {
        return qd_real::_eps;
    }


};
}

#endif /* FEELPP_HAS_QD */

#if defined( FEELPP_HAS_MPFR )
#include <feel/feelcore/mpfr.hpp>

namespace Feel
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
    mpfr_set_default_prec ( __prec );
}

const mp_type mp_eps = mpfr::pow( mp_type(  2 ), -mp_type::GetDefaultPrecision()+1 );

}
#endif // FEELPP_HAS_MPFR



#if !defined(MPI_INT64_T)
#define MPI_INT64_T MPI_LONG_INT
#endif

#if !defined(MPI_INT32_T)
#define MPI_INT32_T MPI_INT
#endif

#if defined(FEELPP_HAS_OPENMP)
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

#endif /* FEELPP_HAS_OPENMP */

#if !defined( DVLOG_IF )

#ifndef NDEBUG
#define DVLOG_IF(verboselevel, condition) VLOG(verboselevel)
#else
#define DVLOG_IF(verboselevel,condition)                                \
    (true || ( !VLOG_IS_ON(verboselevel) && !(condition))) ?            \
    (void) 0 : google::LogMessageVoidify() & LOG(INFO)
#endif // NDEBUG

#endif // DVLOG_IF


#include <feel/feelcore/ptr.hpp>

#endif
