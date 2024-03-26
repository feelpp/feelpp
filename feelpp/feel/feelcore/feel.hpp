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
//! @date 05 Feb 2017
//! @copyright 2017 Feel++ Consortium
//!
# ifndef __cplusplus
# error You must use C++ for Feel
# endif

# ifndef _FEELPP_HH_
# define _FEELPP_HH_

//! @defgroup Feelpp
//! Feel++ classes and methods

//! @defgroup Core
//! @ingroup Feelpp
//! Core classes provided by the library

//! @defgroup Mesh
//! @ingroup Feelpp
//! Mesh classes and algorithms provided by the library

//! @defgroup Discretization
//! @ingroup Feelpp
//! Discretization classes and algorithms provided by the library

//! @defgroup Filters
//! @ingroup Feelpp
//! Filter classes provided by the library

//! @defgroup DSEL-Variational-Formulation
//! @ingroup Feelpp
//! Variational forms provided by the library

//! @defgroup SpaceTime
//! @ingroup Feelpp
//! Time stepping including space provided by the library

//! @defgroup Timing
//! @ingroup Feelpp
//! Timing methods provided by the library

//! @defgroup Traits
//! @ingroup Feelpp
//! Traits provided by the library

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

#include <boost/mp11/integral.hpp>
#include <boost/mp11/utility.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/function.hpp>

#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <boost/assign/list_of.hpp>

#include <boost/math/constants/constants.hpp>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdivision-by-zero"
#pragma clang diagnostic ignored "-Wexpansion-to-defined"
#endif
#include <boost/mpi.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <boost/program_options.hpp>

#include <boost/cstdint.hpp>

#include <boost/format.hpp>
#include <boost/foreach.hpp>

#include <feel/feelcore/hana.hpp>
#include <boost/ref.hpp>

#include <boost/property_tree/ptree.hpp>

#include <filesystem>
#include <fstream>
#include <cmath>
#include <numeric>
#include <string>
#include <limits>
#include <iosfwd>

#if defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning(disable:780)
#endif
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#endif
#if defined(__GNUC__) && !(defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcpp"
#endif
#include <glog/logging.h>
#include <glog/stl_logging.h>
#if defined(__GNUC__) && !(defined(__clang__))
#pragma GCC diagnostic pop
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif

#include <feel/feelconfig.h>
//#include <feel/feelcore/info.hpp>
#include <feel/feelcore/feelmacros.hpp>
#include <feel/feelcore/feelassert.hpp>
#include <feel/feelcore/feeltypes.hpp>
#include <feel/feelcore/feelmath.hpp>

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
namespace fs = std::filesystem;
namespace mpl = boost::mpl;
namespace mp11 = boost::mp11;
namespace po = boost::program_options;
namespace hana=boost::hana;
using namespace boost::hana::literals;
namespace pt =  boost::property_tree;


// bring boost.mpi into Feel realm
namespace mpi=boost::mpi;

namespace constants=boost::math::constants;
//namespace constants = boost::math::constants;
//using boost::math::double_constants;
inline const double pi = constants::pi<double>();
inline const double two_pi = constants::two_pi<double>();

namespace algorithm=boost::algorithm;
using google::WARNING;
using google::ERROR;
using google::INFO;
using google::FATAL;
using boost::format;

using boost::unwrap_ref;

bool filename_is_dot( fs::path const& p );
bool filename_is_dot_dot( fs::path const& p );

//! @cond
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
//! @endcond

} // end namespace Feel

#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Feel
{
//!
//! \namespace
//! @brief Feel++ alias for program_options namespace
//!
namespace po = boost::program_options;

//!
//! @ingroup Core
//! Function that compute the command line option prefix for the option variable map
//!
//! @param prefix prefix
//! @param opt option
//! @param sep separator between sections of the option
//!
std::string
prefixvm( std::string const& prefix,
          std::string const& opt,
          std::string const& sep="." );

//!
//! @ingroup Core
//! @brief trim string to remove special characters
//!
//! trim string removing all leading and trailing spaces and replace
//! all special characters " ;:," inside the block by a _
//! @param s a string to be trimmed
//!
std::string sanitize( std::string const& s );

//!
//! trim a vector of strings removing all leading and trailing spaces and
//! replace all special characters " ;:," inside the block by a _
//! @param s a vector of strings
//!
std::vector<std::string> sanitize( std::vector<std::string> const& s );

//!
//! \namespace
//! @ingroup Core
//!
//! Feel++ namespace alias for boost::posix_time
namespace posix_time = boost::posix_time;
//! Feel++ namespace alias for boost::gregorian
namespace gregorian = boost::gregorian;

//! @cond
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
//! @endcond
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


//! @cond
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
//! @endcond
#endif /* FEELPP_HAS_QD */

#if defined( FEELPP_HAS_MPFR )
#include <feel/feelcore/mpfr.hpp>

//! @cond
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
//! @endcond
#endif // FEELPP_HAS_MPFR

#if defined(FEELPP_DOXYGEN_INVOKED)

//!
//! @ingroup Core
//! Log if condition satisfied when `NDEBUG` is not defined
//!
//! This macro is enabled only if `NDEBUG` is not defined which is to say that
//! the macro will be activited when debugging
//!
//! @code
//! DVLOG_IF( INFO, param < 10 ) << "print only if param < 10";
//! @endcode
//!
#define DVLOG_IF(verboselevel, condition) unspecified

#else

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

# endif // FEELPP_DOXYGEN_INVOKED

/**
 * @brief enable reduce for mpi communication
 */
constexpr inline bool do_reduce = true;
/**
 * @brief disable reduce for mpi communication
 */
constexpr inline bool no_reduce = false;
/**
 * @brief enable communication for mpi communication
 */
constexpr inline bool do_communication = true;
/**
 * @brief disable communication for mpi communication
 */
constexpr inline bool no_communication = false;

/**
 * @brief enable communication for mpi communication
 */
constexpr inline bool parallelEvaluation = true;

/**
 * @brief disable communication for mpi communication
 */
constexpr inline bool sequentialEvaluation = false;

#include <feel/feelcore/ptr.hpp>
#include <feel/feelcore/range.hpp>
#include <feel/feelcore/hashtables.hpp>

/**
 * @brief boolean to enable shared_from_this
 */
inline constexpr bool EnableSharedFromThis = true;
/**
 * @brief boolean to disable shared_from_this
 */
inline constexpr bool DisableSharedFromThis = false;

#endif
