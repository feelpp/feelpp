/** -*- c++ -*- \file cholesky.hpp \brief cholesky decomposition */
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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>


namespace ublas = boost::numeric::ublas;


/** \brief decompose the symmetric positive definit matrix A into product L L^T.
 *
 * template parameter \p MATRIX provides the  type of input matrix
 * template parameter \p TRIA provides the type of lower triangular output matrix
 *
 * \param A square symmetric positive definite input matrix (only the lower triangle is accessed)
 * \param L lower triangular output matrix
 * \return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
 */
template < class MATRIX, class TRIA >
size_t cholesky_decompose( const MATRIX& A, TRIA& L )
{
    using namespace ublas;

    typedef typename MATRIX::value_type T;

    assert( A.size1() == A.size2() );
    assert( A.size1() == L.size1() );
    assert( A.size2() == L.size2() );

    const size_t n = A.size1();

    for ( size_t k=0 ; k < n; k++ )
    {

        double qL_kk = A( k,k ) - inner_prod( project( row( L, k ), range( 0, k ) ),
                                              project( row( L, k ), range( 0, k ) ) );

        if ( qL_kk <= 0 )
        {
            return 1 + k;
        }

        else
        {
            double L_kk = sqrt( qL_kk );
            L( k,k ) = L_kk;

            matrix_column<TRIA> cLk( L, k );
            project( cLk, range( k+1, n ) )
                = ( project( column( A, k ), range( k+1, n ) )
                    - prod( project( L, range( k+1, n ), range( 0, k ) ),
                            project( row( L, k ), range( 0, k ) ) ) ) / L_kk;
        }
    }

    return 0;
}
