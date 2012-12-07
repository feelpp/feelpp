/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-08-17

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
   \file glas.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-08-17
 */
#ifndef __FEELPP_GLAS_HPP
#define __FEELPP_GLAS_HPP 1

#include <string>
#include <sstream>
#include <fstream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/algorithm.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/debug.hpp>
#include <feel/feelcore/traits.hpp>

namespace Feel
{
const double Pi = 3.14159265358979323846264338328;
const double TGV = 1e20;

template <class T>
inline T Min ( const T &a, const T &b )
{
    return a < b ? a : b;
}
template <class T>
inline T Max ( const T &a, const T & b )
{
    return a > b ? a : b;
}
template <class T>
inline T Abs ( const T &a )
{
    return a < 0 ? -a : a;
}

template <class T>
inline void Exchange ( T& a, T& b )
{
    T c = a;
    a = b;
    b = c;
}
template <class T>
inline T Max ( const T &a, const T & b, const T & c )
{
    return Max( Max( a, b ), c );
}
template <class T>
inline T Min ( const T &a, const T & b, const T & c )
{
    return Min( Min( a, b ), c );
}

namespace ublas = boost::numeric::ublas;



/* reduce operations */
using ublas::sum;
using ublas::norm_1;
using ublas::norm_2;
using ublas::norm_inf;
using ublas::index_norm_inf;

/* binary operations */
using ublas::outer_prod;
using ublas::inner_prod;

struct norm_inf_adaptor
{
    template<typename E>
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
    //typedef ublas::vector<T, ublas::bounded_array<T, S> >  type;
    typedef ublas::vector<T>  type;
};
#if 0
/*!
  \typedef Vector node_type
  type for nodes. Typically \c node_type will be of size 1,2,...,3
*/
typedef node<double, 3>::type node_type;

/*!
  gradient type
*/
typedef node<double, 3>::type  gradient_node_type;
#endif
/*!
  hessian type
*/
typedef ublas::symmetric_matrix<double, ublas::lower, ublas::row_major, ublas::bounded_array<double, 9> >  hessian_node_type;


typedef ublas::matrix<double, ublas::column_major, ublas::bounded_array<double, 9> >  lapack_matrix_type;
typedef ublas::symmetric_adaptor<lapack_matrix_type, ublas::lower> symmetric_matrix_type;

typedef lapack_matrix_type transformation_matrix_type;

/**
 * matrix node definition
 */
template <typename T = double, uint16_type S = 256>
struct matrix_node
{
    //typedef ublas::matrix<T, ublas::column_major, ublas::bounded_array<T, S> >  type;
    typedef ublas::matrix<T, ublas::column_major>  type;
};

#if 0
typedef matrix_node<double>::type matrix_node_type;
#endif

inline
DebugStream&
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
inline
NdebugStream&
operator<<( NdebugStream& os, node<real64_type>::type const& /*n*/ )
{
    return os;
}


#if defined( FEELPP_HAS_QD_QD_H )
inline
DebugStream&
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
inline
NdebugStream&
operator<<( NdebugStream& __os, node<dd_real>::type const& __n )
{
    return __os;
}

inline
DebugStream&
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
inline
NdebugStream&
operator<<( NdebugStream& __os, node<qd_real>::type const& __n )
{
    return __os;
}


#endif /* FEELPP_HAS_QD_QD_H */

template<typename T>
inline
DebugStream&
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
template<typename T>
inline
NdebugStream&
operator<<( NdebugStream& __os, ublas::vector<T> const& /*__n*/ )
{
    return __os;
}
template<typename T, typename Orient>
inline
DebugStream&
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
template<typename T,typename Orient>
inline
NdebugStream&
operator<<( NdebugStream& __os, ublas::matrix<T,Orient> const& /*__n*/ )
{
    return __os;
}


//
// sparse matrices
//
#if defined( FEELPP_SIZET_SAME_AS_UINT )
typedef ublas::compressed_matrix<double,
        ublas::row_major > csr_matrix_type;
typedef ublas::compressed_matrix<double,
        ublas::column_major > csc_matrix_type;
#else
typedef ublas::compressed_matrix<double,
        ublas::row_major, 0,
        ublas::unbounded_array<int>,
        ublas::unbounded_array<double> > csr_matrix_type;
typedef ublas::compressed_matrix<double,
        ublas::column_major, 0,
        ublas::unbounded_array<int>,
        ublas::unbounded_array<double> > csc_matrix_type;
#endif

/**
 * Dump vector to file in Matlab format and spy
 *
 */
template<typename MatrixType>
void spy( MatrixType const& __m, std::string const &filename )
{
    std::string name = filename;
    std::string separator = " , ";

    // check on the file name
    int i = filename.find( "." );

    if ( i <= 0 )
        name = filename + ".m";

    else
    {
        if ( ( unsigned int ) i != filename.size() - 2 ||
                filename[ i + 1 ] != 'm' )
        {
            std::cerr << "Wrong file name ";
            name = filename + ".m";
        }
    }

    std::ofstream file_out( name.c_str() );

    FEELPP_ASSERT( file_out )( filename ).error( "[Feel::spy] ERROR: File cannot be opened for writing." );

    file_out << "S = [ ";

    for ( typename MatrixType::const_iterator1 i1=__m.begin1();
            i1!=__m.end1(); ++i1 )
    {
        for ( typename MatrixType::const_iterator2 i2=i1.begin();
                i2!=i1.end(); ++i2 )
            file_out << i2.index1() + 1 << separator
                     << i2.index2() + 1 << separator
                     << *i2  << std::endl;
    }

    file_out << "];" << std::endl;
    file_out << "I=S(:,1); J=S(:,2); S=S(:,3);" << std::endl;
    file_out << "A=sparse(I,J,S); spy(A);" << std::endl;
}

namespace glas
{
namespace ublas = boost::numeric::ublas;

template<typename T, typename Orien>
inline
ublas::matrix<T,Orien>
average( ublas::matrix<T, Orien> const& m )
{
    ublas::matrix<T, Orien> v( m.size1(), 1 );
    ublas::scalar_vector<T> avg( m.size2(), T( 1 ) );
    T n_val = int( m.size2() );

    for ( size_type i = 0; i < v.size1(); ++i )
        v( i, 0 ) = ublas::inner_prod( ublas::row( m, i ), avg )/n_val;

    return v;
}
template<typename T>
inline
void
clean( T& t,
       typename T::value_type const& treshold = type_traits<typename T::value_type>::epsilon(),
       typename T::value_type const& new_value = typename T::value_type( 0.0 ) )
{
    std::for_each( t.data().begin(),
                   t.data().end(),
                   lambda::if_then( lambda::_1  < lambda::constant( treshold ) && lambda::_1  > -lambda::constant( treshold ),
                                    lambda::_1 = lambda::constant( new_value ) ) );
}

template<typename T>
inline
ublas::vector<T>
linspace( T const& __a, T const& __b, size_type __N, int interior = 0 )
{
    size_type N = __N-2*interior;
    ublas::vector<T> v( N );
    T h = ( __b-__a )/T( __N-1 );
    T a = __a+T( interior )*h;
    //T b = __b-T(interior)*h;
#if 0
    size_type i = 0;
    std::for_each( v.begin(),
                   v.end(),
                   ( lambda::_1 = a+lambda::var( i )*h, ++lambda::var( i ) ) );
#else

    for ( size_type i = 0; i < N; ++i )
        v[i]=a+T( i )*h;

#endif
    return v;
}

template<typename T>
inline
void
randomize( T& t )
{
    typedef typename T::value_type value_type;
    typedef boost::minstd_rand base_generator_type;
    base_generator_type generator( 42u );
    boost::uniform_real<value_type> uni_dist( 0.0,1.0 );
    typedef boost::variate_generator<base_generator_type&, boost::uniform_real<value_type> >  rand_type;
    rand_type uni( generator, uni_dist );

    std::for_each( t.data().begin(),
                   t.data().end(),
                   lambda::_1 = lambda::bind<value_type>( uni ) );

}

//
//
//
} // namespace glas
} // namespace Feel

#endif /* __FEELPP_GLAS_HPP */
