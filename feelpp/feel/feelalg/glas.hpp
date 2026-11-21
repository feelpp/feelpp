/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2005, 2006 EPFL
    SPDX-FileCopyrightText: 2011-2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
   \file glas.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-17
 */
#ifndef __FEELPP_GLAS_HPP
#define __FEELPP_GLAS_HPP 1

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <Eigen/src/Core/util/ReenableStupidWarnings.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/debug.hpp>
#include <feel/feelcore/traits.hpp>

#include <feel/feelalg/eigenmap.hpp>
#include <feel/feelalg/fmt.hpp>

namespace Feel
{
const double Pi = 3.14159265358979323846264338328;
const double TGV = 1e20;

/**
 * @return true is the matrix has no Inf, false otherwise
 */
template <typename Derived>
inline bool is_finite( const Eigen::MatrixBase<Derived>& x )
{
    return ( ( x - x ).array() == ( x - x ).array() ).all();
}

/**
 * @return true is matrix has no NaN, false otherwise
 */
template <typename Derived>
inline bool is_nan( const Eigen::MatrixBase<Derived>& x )
{
    return ( ( x.array() == x.array() ) ).all();
}

template <class T>
inline T Min( const T& a, const T& b )
{
    return a < b ? a : b;
}
template <class T>
inline T Max( const T& a, const T& b )
{
    return a > b ? a : b;
}
template <class T>
inline T Abs( const T& a )
{
    return a < 0 ? -a : a;
}

template <class T>
inline void Exchange( T& a, T& b )
{
    T c = a;
    a = b;
    b = c;
}
template <class T>
inline T Max( const T& a, const T& b, const T& c )
{
    return Max( Max( a, b ), c );
}
template <class T>
inline T Min( const T& a, const T& b, const T& c )
{
    return Min( Min( a, b ), c );
}

namespace ublas = boost::numeric::ublas;

/* reduce operations */
using ublas::index_norm_inf;
using ublas::norm_1;
using ublas::norm_2;
using ublas::norm_inf;
using ublas::sum;

/* binary operations */
using ublas::inner_prod;
using ublas::outer_prod;

struct norm_inf_adaptor
{
    template <typename E>
    Real operator()( E const& __v ) const
    {
        return norm_inf( __v );
    }
};

/**
 * node definition
 */
template <typename T = double, uint16_type S = 3>
struct node
{
    // typedef ublas::vector<T, ublas::bounded_array<T, S> >  type;
    typedef ublas::vector<T> type;
};

template <typename T = double, uint16_type S = 3>
using node_t = typename node<T, S>::type;

/*!
  hessian type
*/
typedef ublas::symmetric_matrix<double, ublas::lower, ublas::row_major, ublas::bounded_array<double, 9>> hessian_node_type;

typedef ublas::matrix<double, ublas::column_major, ublas::bounded_array<double, 9>> lapack_matrix_type;
typedef ublas::symmetric_adaptor<lapack_matrix_type, ublas::lower> symmetric_matrix_type;

typedef lapack_matrix_type transformation_matrix_type;

/**
 * matrix node definition
 */
template <typename T = double, uint16_type S = 256>
struct matrix_node
{
    // typedef ublas::matrix<T, ublas::column_major, ublas::bounded_array<T, S> >  type;
    typedef ublas::matrix<T, ublas::column_major> type;
};

template <typename T = double>
using matrix_node_t = ublas::matrix<T, ublas::column_major>;

//! Eigen type to map matrix_type (row major)
template <typename T = double>
using em_matrix_row_type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

//! Eigen type to map matrix_node<>::type
template <typename T = double>
using em_matrix_col_type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;

//! Eigen type to map matrix_node<>::type
template <typename T = double>
using em_cmatrix_col_type = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;

/**
 * @brief Compatibility spelling for a mutable zero-copy Eigen matrix view.
 * @tparam T scalar type
 * @tparam Layout uBLAS matrix layout
 * @tparam Storage uBLAS contiguous storage
 * @param m owning uBLAS matrix
 * @return mutable Eigen view
 * @deprecated Prefer the explicit `eigenMap()` spelling in new code.
 */
template <typename T, typename Layout, typename Storage>
auto emap( ublas::matrix<T, Layout, Storage>& m ) noexcept
{
    return eigenMap( m );
}
/**
 * @brief Compatibility spelling for a read-only zero-copy Eigen matrix view.
 * @tparam T scalar type
 * @tparam Layout uBLAS matrix layout
 * @tparam Storage uBLAS contiguous storage
 * @param m owning uBLAS matrix
 * @return read-only Eigen view
 * @deprecated Prefer the explicit `eigenMap()` spelling in new code.
 */
template <typename T, typename Layout, typename Storage>
auto emap( ublas::matrix<T, Layout, Storage> const& m ) noexcept
{
    return eigenMap( m );
}

//! Eigen type to map node_type
template <typename T = double>
using em_node_type = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>;

/**
 * @brief Compatibility spelling for a mutable zero-copy Eigen vector view.
 * @tparam T scalar type
 * @tparam Storage uBLAS contiguous storage
 * @param m owning uBLAS vector
 * @return mutable Eigen column-vector view
 * @deprecated Prefer the explicit `eigenMap()` spelling in new code.
 */
template <typename T, typename Storage>
auto emap( ublas::vector<T, Storage>& m ) noexcept
{
    return eigenMap( m );
}
/**
 * @brief Compatibility spelling for a read-only zero-copy Eigen vector view.
 * @tparam T scalar type
 * @tparam Storage uBLAS contiguous storage
 * @param m owning uBLAS vector
 * @return read-only Eigen column-vector view
 * @deprecated Prefer the explicit `eigenMap()` spelling in new code.
 */
template <typename T, typename Storage>
auto emap( ublas::vector<T, Storage> const& m ) noexcept
{
    return eigenMap( m );
}

//!
//! Tensor map to a ublas_type
//!
template <typename value_type = double>
using tensor_eigen_ublas_type = Eigen::TensorMap<Eigen::Tensor<const value_type, 2>>;

using d_tensor_eigen_ublas_type = Eigen::TensorMap<Eigen::Tensor<const double, 2>>;

template <int Dim, typename T = double>
using eigen_vector_type = Eigen::Matrix<T, Dim, 1>;

template <typename T = double>
using eigen_vector_x_col_type = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <int Dim, typename T = double>
using eigen_vector_x_row_type = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template <int Dim, typename T = double>
using vector_eigen_vector_type = std::vector<eigen_vector_type<Dim, T>, Eigen::aligned_allocator<eigen_vector_type<Dim, T>>>;

template <int N = Eigen::Dynamic, int P = Eigen::Dynamic, typename T = double>
using eigen_matrix_type = Eigen::Matrix<T, N, P>;

template <typename T = double>
using eigen_matrix_xx_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T = double>
using eigen_matrix_xx_col_type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

template <int N, int P, typename T = double>
using vector_eigen_matrix_type = std::vector<eigen_matrix_type<N, P, T>, Eigen::aligned_allocator<eigen_matrix_type<N, P, T>>>;

template <typename Key, int N, int P, typename T = double>
using allocator_matrix_t = Eigen::aligned_allocator<std::pair<const Key, eigen_matrix_type<N, P, T>>>;

template <typename Key, int N, int P, typename T = double>
using unordered_map_eigen_matrix_type = std::unordered_map<Key, eigen_matrix_type<N, P, T>, boost::hash<Key>, std::equal_to<Key>, allocator_matrix_t<Key, N, P, T>>;

template <int M, int N, typename value_type = double>
using tensor2_fixed_size_t = Eigen::TensorFixedSize<value_type, Eigen::Sizes<M, N>>;

template <int M, int N, int P, typename value_type = double>
using tensor3_fixed_size_t = Eigen::TensorFixedSize<value_type, Eigen::Sizes<M, N, P>>;

template <int M, int N, typename value_type = double>
using vector_tensor2_fixed_size_t = std::vector<tensor2_fixed_size_t<M, N, value_type>>;

template <int M, int N, int P, typename value_type = double>
using vector_tensor3_fixed_size_t = std::vector<tensor3_fixed_size_t<M, N, P, value_type>>;
//!
//! Tensor Maps
//!
template <int M, int N, typename value_type = double>
using tensor_map_fixed_size_matrix_t = Eigen::TensorMap<Eigen::TensorFixedSize<const value_type, Eigen::Sizes<M, N>>>;

template <int N, typename value_type = double>
using tensor_map_t = Eigen::TensorMap<Eigen::Tensor<const value_type, N>>;

template <typename value_type = double>
using tensor2_map_t = Eigen::TensorMap<Eigen::Tensor<const value_type, 2>>;

template <typename value_type = double>
using tensor1_map_t = Eigen::TensorMap<Eigen::Tensor<const value_type, 1>>;

template <int M, int N, typename value_type = double>
using em_fixed_size_matrix_t = Eigen::Map<Eigen::Matrix<value_type, M, N>>;
template <int M, int N, typename value_type = double>
using em_fixed_size_cmatrix_t = Eigen::Map<const Eigen::Matrix<value_type, M, N>>;

//!
//! dimension array for Eigen::Tensor
//!
using dimpair_t = Eigen::IndexPair<int>;

inline FEELPP_DEPRECATED DebugStream&
operator<<( DebugStream& __os, node<real64_type>::type const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}
inline FEELPP_DEPRECATED NdebugStream&
operator<<( NdebugStream& os, node<real64_type>::type const& /*n*/ )
{
    return os;
}

#if defined( FEELPP_HAS_QD_QD_H )
inline DebugStream& FEELPP_DEPRECATED
operator<<( DebugStream& __os, node<dd_real>::type const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}
inline FEELPP_DEPRECATED NdebugStream&
operator<<( NdebugStream& __os, node<dd_real>::type const& __n )
{
    return __os;
}

inline DebugStream& FEELPP_DEPRECATED
operator<<( DebugStream& __os, node<qd_real>::type const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}
inline FEELPP_DEPRECATED NdebugStream&
operator<<( NdebugStream& __os, node<qd_real>::type const& __n )
{
    return __os;
}

#endif /* FEELPP_HAS_QD_QD_H */

template <typename T>
inline FEELPP_DEPRECATED DebugStream&
operator<<( DebugStream& __os, ublas::vector<T> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}
template <typename T>
inline FEELPP_DEPRECATED NdebugStream&
operator<<( NdebugStream& __os, ublas::vector<T> const& /*__n*/ )
{
    return __os;
}
template <typename T, typename Orient>
inline FEELPP_DEPRECATED DebugStream&
operator<<( DebugStream& __os, ublas::matrix<T, Orient> const& __n )
{
    if ( __os.doPrint() )
    {
        std::ostringstream __str;

        __str << __n;

        __os << __str.str() << "\n";
    }

    return __os;
}
template <typename T, typename Orient>
inline FEELPP_DEPRECATED NdebugStream&
operator<<( NdebugStream& __os, ublas::matrix<T, Orient> const& /*__n*/ )
{
    return __os;
}

//
// sparse matrices
//
#if defined( FEELPP_SIZET_SAME_AS_UINT )
typedef ublas::compressed_matrix<double,
                                 ublas::row_major>
    csr_matrix_type;
typedef ublas::compressed_matrix<double,
                                 ublas::column_major>
    csc_matrix_type;
#else
typedef ublas::compressed_matrix<double,
                                 ublas::row_major, 0,
                                 ublas::unbounded_array<int>,
                                 ublas::unbounded_array<double>>
    csr_matrix_type;
typedef ublas::compressed_matrix<double,
                                 ublas::column_major, 0,
                                 ublas::unbounded_array<int>,
                                 ublas::unbounded_array<double>>
    csc_matrix_type;
#endif

/**
 * Dump vector to file in Matlab format and spy
 *
 */
template <typename MatrixType>
void spy( MatrixType const& __m, std::string const& filename )
{
    std::string name = filename;
    std::string separator = " , ";

    // check on the file name
    int i = filename.find( "." );

    if ( i <= 0 )
        name = filename + ".m";

    else
    {
        if ( (unsigned int)i != filename.size() - 2 ||
             filename[i + 1] != 'm' )
        {
            std::cerr << "Wrong file name ";
            name = filename + ".m";
        }
    }

    std::ofstream file_out( name.c_str() );

    FEELPP_ASSERT( file_out )
    ( filename ).error( "[Feel::spy] ERROR: File cannot be opened for writing." );

    file_out << "S = [ ";

    for ( typename MatrixType::const_iterator1 i1 = __m.begin1();
          i1 != __m.end1(); ++i1 )
    {
        for ( typename MatrixType::const_iterator2 i2 = i1.begin();
              i2 != i1.end(); ++i2 )
            file_out << i2.index1() + 1 << separator
                     << i2.index2() + 1 << separator
                     << *i2 << std::endl;
    }

    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3);" << std::endl;
    file_out << "A=sparse(I,J,S); spy(A);" << std::endl;
}

namespace glas
{
namespace ublas = boost::numeric::ublas;

template <typename T, typename Orien>
inline ublas::matrix<T, Orien>
average( ublas::matrix<T, Orien> const& m )
{
    ublas::matrix<T, Orien> v( m.size1(), 1 );
    ublas::scalar_vector<T> avg( m.size2(), T( 1 ) );
    T n_val = int( m.size2() );

    for ( size_type i = 0; i < v.size1(); ++i )
        v( i, 0 ) = ublas::inner_prod( ublas::row( m, i ), avg ) / n_val;

    return v;
}
template <typename T>
inline void
clean( T& t,
       typename T::value_type const& treshold = type_traits<typename T::value_type>::epsilon(),
       typename T::value_type const& new_value = typename T::value_type( 0.0 ) )
{
    std::for_each( t.data().begin(),
                   t.data().end(),
                   [new_value, treshold]( auto& t )
                   {
                       if ( std::abs( t ) < treshold )
                           t = new_value;
                   } );
}

template <typename T>
inline ublas::vector<T>
linspace( T const& __a, T const& __b, size_type __N, int interior = 0 )
{
    size_type N = __N - 2 * interior;
    ublas::vector<T> v( N );
    T h = ( __b - __a ) / T( __N - 1 );
    T a = __a + T( interior ) * h;
    // T b = __b-T(interior)*h;
#if 0
    size_type i = 0;
    std::for_each( v.begin(),
                   v.end(),
                   ( lambda::_1 = a+lambda::var( i )*h, ++lambda::var( i ) ) );
#else

    for ( size_type i = 0; i < N; ++i )
        v[i] = a + T( i ) * h;

#endif
    return v;
}

template <typename T>
inline void
randomize( T& t )
{
    typedef typename T::value_type value_type;
    typedef boost::minstd_rand base_generator_type;
    base_generator_type generator( 42u );
    boost::uniform_real<value_type> uni_dist( 0.0, 1.0 );
    typedef boost::variate_generator<base_generator_type&, boost::uniform_real<value_type>> rand_type;
    rand_type uni( generator, uni_dist );

    std::for_each( t.data().begin(),
                   t.data().end(),
                   [&uni]( auto& t )
                   { t = uni; } );
}

//
//
//
} // namespace glas
} // namespace Feel

#endif /* __FEELPP_GLAS_HPP */
