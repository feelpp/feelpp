/** -*- c++ -*- \file cholesky_test.cpp \brief test cholesky decomposition */
/*
 -   begin                : 2005-08-24
 -   copyright            : (C) 2005 by Gunter Winkler, Konstantin Kutzkow
 -   email                : guwi17@gmx.de

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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <cassert>
#include <limits>

#include <boost/timer.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/io.hpp>

#include "cholesky.hpp"

namespace ublas = boost::numeric::ublas;

/** \brief make a immutable triangular adaptor from a matrix
 *
 * \usage:
<code>
 A = triangular< lower >(B);
 A = triangular(B, lower());
</code>
 */
template < class TYPE, class MATRIX >
ublas::triangular_adaptor<const MATRIX, TYPE>
triangular( const MATRIX & A, const TYPE& uplo = TYPE() )
{
    return ublas::triangular_adaptor<const MATRIX, TYPE>( A );
}

/** \brief make a immutable banded adaptor from a matrix
 *
 * \usage:
<code>
 A = banded(B, lower, upper);
</code>
 */
template < class MATRIX >
ublas::banded_adaptor<const MATRIX>
banded( const MATRIX & A, const size_t lower, const size_t upper )
{
    return ublas::banded_adaptor<const MATRIX>( A, lower, upper );
}

/** \brief fill lower triangular matrix L
 */
template < class MATRIX >
void fill_symm( MATRIX & L, const size_t bands = std::numeric_limits<size_t>::max() )
{
    typedef typename MATRIX::size_type size_type;

    assert( L.size1() == L.size2() );

    size_type size = L.size1();

    for ( size_type i=0; i<size; i++ )
    {
        for ( size_type j = ( ( i>bands )?( i-bands ):0 ); j<i; j++ )
        {
            L( i,j ) = 1 + ( ( double ) i )*j + ( ( double ) j )*j;
        }

        L( i,i ) = 20*size;
    }

    return;
}

int main( int argc, char * argv[] )
{
    size_t size = 10;

    if ( argc > 1 )
        size = ::atoi ( argv [1] );

    boost::timer  t1;
    double pr, de;

    typedef double DBL;
    {
        // use dense matrix
        ublas::matrix<DBL> A ( size, size );
        ublas::matrix<DBL> T ( size, size );
        ublas::matrix<DBL> L ( size, size );

        A = ublas::zero_matrix<DBL>( size, size );

        fill_symm( T );
        t1.restart();
        A = ublas::prod( T, trans( T ) );
        pr = t1.elapsed();
        t1.restart();
        size_t res = cholesky_decompose( A, L );
        de = t1.elapsed();

        std::cout << res << ": "
                  << ublas::norm_inf( L-T )
                  << " (deco: " << de << " sec)"
                  << " (prod: " << pr << " sec)"
                  << " " << size
                  << std::endl;
    }

    {
        // use dense triangular matrices
        ublas::triangular_matrix<DBL, ublas::lower> A ( size, size );
        ublas::triangular_matrix<DBL, ublas::lower> T ( size, size );
        ublas::triangular_matrix<DBL, ublas::lower> L ( size, size );

        A = ublas::zero_matrix<DBL> ( size, size ) ;
        A = triangular<ublas::lower>( ublas::zero_matrix<DBL> ( size, size ) );
        A = triangular( ublas::zero_matrix<DBL> ( size, size ), ublas::lower() );

        fill_symm( T );
        t1.restart();
        A = triangular<ublas::lower>( ublas::prod( T, trans( T ) ) );
        pr = t1.elapsed();
        t1.restart();
        size_t res = cholesky_decompose( A, L );
        de = t1.elapsed();

        std::cout << res << ": "
                  << ublas::norm_inf( L-T )
                  << " (deco: " << de << " sec)"
                  << " (prod: " << pr << " sec)"
                  << " " << size
                  << std::endl;
    }

    {
        // use banded matrices matrix
        typedef ublas::banded_matrix<DBL> MAT;

        size_t bands = std::min<size_t>( size, 50 );
        MAT A ( size, size, bands, 0 );
        MAT T ( size, size, bands, 0 );
        MAT L ( size, size, bands, 0 );

        A = ublas::zero_matrix<DBL> ( size, size ) ;
        A = banded( ublas::zero_matrix<DBL> ( size, size ), bands, 0 );

        fill_symm( T, bands );
        t1.restart();
        A = banded( ublas::prod( T, trans( T ) ), bands, 0 );
        pr = t1.elapsed();
        t1.restart();
        size_t res = cholesky_decompose( A, L );
        de = t1.elapsed();

        std::cout << res << ": "
                  << ublas::norm_inf( L-T )
                  << " (deco: " << de << " sec)"
                  << " (prod: " << pr << " sec)"
                  << " " << size
                  << std::endl;
    }

    return EXIT_SUCCESS;
}


/****************

$ g++ -I $HOME/include -o cholesky_test cholesky_test.cpp -DNDEBUG -O2

$ for i in 100 200 400 800 1600; do ./cholesky_test $i ; done
0: 0 (deco: 0 sec) (prod: 0.02 sec) 100
0: 0 (deco: 0.01 sec) (prod: 0.01 sec) 100
0: 0 (deco: 0.01 sec) (prod: 0.01 sec) 100
0: 0 (deco: 0.02 sec) (prod: 0.13 sec) 200
0: 0 (deco: 0.05 sec) (prod: 0.06 sec) 200
0: 0 (deco: 0.06 sec) (prod: 0.06 sec) 200
0: 0 (deco: 0.18 sec) (prod: 1.03 sec) 400
0: 0 (deco: 0.43 sec) (prod: 0.42 sec) 400
0: 0 (deco: 0.47 sec) (prod: 0.5 sec) 400
0: 0 (deco: 1.55 sec) (prod: 8.16 sec) 800
0: 0 (deco: 3.41 sec) (prod: 3.2 sec) 800
0: 0 (deco: 4.6 sec) (prod: 4.87 sec) 800
0: 0 (deco: 12.14 sec) (prod: 65.15 sec) 1600
0: 0 (deco: 26.74 sec) (prod: 25.47 sec) 1600
0: 0 (deco: 70.29 sec) (prod: 75.68 sec) 1600

(note: dense + tria + full banded)

 ****************/
