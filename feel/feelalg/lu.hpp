/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-08

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
   \file lu.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-08
 */
#ifndef __DenseLU_H
#define __DenseLU_H 1

#include <boost/mpl/int.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/glas.hpp>

namespace Feel
{
namespace mpl = boost::mpl;

namespace details
{
template<typename Matrix>
inline typename Matrix::value_type
det( Matrix const& M, mpl::int_<1> )
{
    return M( 0, 0 );
}
template<typename Matrix>
inline typename Matrix::value_type
det( Matrix const& M, mpl::int_<2> )
{
    return  M( 0, 0 )*M( 1, 1 )-M( 0, 1 )*M( 1, 0 );
}
template<typename Matrix>
inline typename Matrix::value_type
det( Matrix const& M, mpl::int_<3> )
{
    return ( M( 0, 0 )*M( 1, 1 )*M( 2, 2 )-
             M( 0, 0 )*M( 1, 2 )*M( 2, 1 )-
             M( 1, 0 )*M( 0, 1 )*M( 2, 2 )+
             M( 1, 0 )*M( 0, 2 )*M( 2, 1 )+
             M( 2, 0 )*M( 0, 1 )*M( 1, 2 )-
             M( 2, 0 )*M( 0, 2 )*M( 1, 1 ) );
}

template<typename Matrix>
inline typename Matrix::value_type
det( Matrix const& M, mpl::int_<4> )
{
    return ( M( 0, 0 ) * M( 1, 1 ) * M( 2, 2 ) * M( 3, 3 ) -
             M( 0, 0 ) * M( 1, 1 ) * M( 2, 3 ) * M( 3, 2 ) -
             M( 0, 0 ) * M( 2, 1 ) * M( 1, 2 ) * M( 3, 3 ) +
             M( 0, 0 ) * M( 2, 1 ) * M( 1, 3 ) * M( 3, 2 ) +
             M( 0, 0 ) * M( 3, 1 ) * M( 1, 2 ) * M( 2, 3 ) -
             M( 0, 0 ) * M( 3, 1 ) * M( 1, 3 ) * M( 2, 2 ) -
             M( 1, 0 ) * M( 0, 1 ) * M( 2, 2 ) * M( 3, 3 ) +
             M( 1, 0 ) * M( 0, 1 ) * M( 2, 3 ) * M( 3, 2 ) +
             M( 1, 0 ) * M( 2, 1 ) * M( 0, 2 ) * M( 3, 3 ) -
             M( 1, 0 ) * M( 2, 1 ) * M( 0, 3 ) * M( 3, 2 ) -
             M( 1, 0 ) * M( 3, 1 ) * M( 0, 2 ) * M( 2, 3 ) +
             M( 1, 0 ) * M( 3, 1 ) * M( 0, 3 ) * M( 2, 2 ) +
             M( 2, 0 ) * M( 0, 1 ) * M( 1, 2 ) * M( 3, 3 ) -
             M( 2, 0 ) * M( 0, 1 ) * M( 1, 3 ) * M( 3, 2 ) -
             M( 2, 0 ) * M( 1, 1 ) * M( 0, 2 ) * M( 3, 3 ) +
             M( 2, 0 ) * M( 1, 1 ) * M( 0, 3 ) * M( 3, 2 ) +
             M( 2, 0 ) * M( 3, 1 ) * M( 0, 2 ) * M( 1, 3 ) -
             M( 2, 0 ) * M( 3, 1 ) * M( 0, 3 ) * M( 1, 2 ) -
             M( 3, 0 ) * M( 0, 1 ) * M( 1, 2 ) * M( 2, 3 ) +
             M( 3, 0 ) * M( 0, 1 ) * M( 1, 3 ) * M( 2, 2 ) +
             M( 3, 0 ) * M( 1, 1 ) * M( 0, 2 ) * M( 2, 3 ) -
             M( 3, 0 ) * M( 1, 1 ) * M( 0, 3 ) * M( 2, 2 ) -
             M( 3, 0 ) * M( 2, 1 ) * M( 0, 2 ) * M( 1, 3 ) +
             M( 3, 0 ) * M( 2, 1 ) * M( 0, 3 ) * M( 1, 2 ) );

}


template<typename Matrix>
inline void
inverse( Matrix const& M, Matrix& Minv, mpl::int_<1> )
{
    Minv( 0, 0 ) = 1.0/M( 0, 0 );
}
template<typename Matrix>
inline void
inverse( Matrix const& M, Matrix& Minv, mpl::int_<2> )
{
    typename Matrix::value_type J = M( 0, 0 )*M( 1, 1 )-M( 0, 1 )*M( 1, 0 );
    Minv( 0, 0 ) = M( 1, 1 )/J;
    Minv( 1, 0 ) = -M( 1, 0 )/J;
    Minv( 0, 1 ) = -M( 0, 1 )/J;
    Minv( 1, 1 ) = M( 0, 0 )/J;
}
template<typename Matrix>
inline void
inverse( Matrix const& M, Matrix& Minv, mpl::int_<3> )
{
    typedef typename Matrix::value_type value_type;

    value_type t4 = M( 0, 0 )*M( 1, 1 );
    value_type t6 = M( 0, 0 )*M( 1, 2 );
    value_type t8 = M( 0, 1 )*M( 1, 0 );
    value_type t10 = M( 0, 2 )*M( 1, 0 );

    value_type t12 = M( 0, 1 )*M( 2, 0 );
    value_type t14 = M( 0, 2 )*M( 2, 0 );
    value_type t17 = 1.0/( M( 2, 2 )*t4-t6*M( 2, 1 )-t8*M( 2, 2 )+t10*M( 2, 1 )+t12*M( 1, 2 )-t14*M( 1, 1 ) );

    Minv( 0, 0 ) = ( M( 1, 1 )*M( 2, 2 )-M( 1, 2 )*M( 2, 1 ) )*t17;
    Minv( 0, 1 ) = -( M( 0, 1 )*M( 2, 2 )-M( 0, 2 )*M( 2, 1 ) )*t17;
    Minv( 0, 2 ) = -( -M( 0, 1 )*M( 1, 2 )+M( 0, 2 )*M( 1, 1 ) )*t17;
    Minv( 1, 0 ) = -( M( 1, 0 )*M( 2, 2 )-M( 1, 2 )*M( 2, 0 ) )*t17;
    Minv( 1, 1 ) = ( M( 0, 0 )*M( 2, 2 )-t14 )*t17;
    Minv( 1, 2 ) = -( t6-t10 )*t17;
    Minv( 2, 0 ) = -( -M( 1, 0 )*M( 2, 1 )+M( 1, 1 )*M( 2, 0 ) )*t17;
    Minv( 2, 1 ) = -( M( 0, 0 )*M( 2, 1 )-t12 )*t17;
    Minv( 2, 2 ) = ( t4-t8 )*t17;
}

template<typename Matrix>
inline void
inverse( Matrix const& __restrict__ M, Matrix& __restrict__ Minv, typename Matrix::value_type const& J, mpl::int_<1> )
{
    Minv( 0, 0 ) = typename Matrix::value_type( 1.0 )/J;
}
template<typename Matrix>
inline void
inverse( Matrix const& __restrict__ M, Matrix& __restrict__ Minv, typename Matrix::value_type const& J, mpl::int_<2> )
{
    Minv( 0, 0 ) = M( 1, 1 )/J;
    Minv( 1, 0 ) = -M( 1, 0 )/J;
    Minv( 0, 1 ) = -M( 0, 1 )/J;
    Minv( 1, 1 ) = M( 0, 0 )/J;
}
template<typename Matrix>
inline void
inverse( Matrix const& __restrict__ M, Matrix& __restrict__ Minv, typename Matrix::value_type const& J, mpl::int_<3> )
{
    typedef typename Matrix::value_type value_type;

    value_type t4 = M( 0, 0 )*M( 1, 1 );
    value_type t6 = M( 0, 0 )*M( 1, 2 );
    value_type t8 = M( 0, 1 )*M( 1, 0 );
    value_type t10 = M( 0, 2 )*M( 1, 0 );

    value_type t12 = M( 0, 1 )*M( 2, 0 );
    value_type t14 = M( 0, 2 )*M( 2, 0 );

    Minv( 0, 0 ) = ( M( 1, 1 )*M( 2, 2 )-M( 1, 2 )*M( 2, 1 ) )/J;
    Minv( 0, 1 ) = -( M( 0, 1 )*M( 2, 2 )-M( 0, 2 )*M( 2, 1 ) )/J;
    Minv( 0, 2 ) = -( -M( 0, 1 )*M( 1, 2 )+M( 0, 2 )*M( 1, 1 ) )/J;
    Minv( 1, 0 ) = -( M( 1, 0 )*M( 2, 2 )-M( 1, 2 )*M( 2, 0 ) )/J;
    Minv( 1, 1 ) = ( M( 0, 0 )*M( 2, 2 )-t14 )/J;
    Minv( 1, 2 ) = -( t6-t10 )/J;
    Minv( 2, 0 ) = -( -M( 1, 0 )*M( 2, 1 )+M( 1, 1 )*M( 2, 0 ) )/J;
    Minv( 2, 1 ) = -( M( 0, 0 )*M( 2, 1 )-t12 )/J;
    Minv( 2, 2 ) = ( t4-t8 )/J;
}


} // details

template<int Dim, typename Matrix>
inline typename Matrix::value_type
det( Matrix const& M )
{
    return details::det( M, mpl::int_<Dim>() );
}

template<int Dim, typename Matrix>
inline void
inverse( Matrix const& M, Matrix& Minv )
{
    details::inverse( M, Minv, mpl::int_<Dim>() );
}

template<int Dim, typename Matrix>
inline void
inverse( Matrix const& __restrict__ M, Matrix& __restrict__ Minv, typename Matrix::value_type const& J )
{
    details::inverse( M, Minv, J, mpl::int_<Dim>() );
}


/** LU Decomposition.
    <P>
    For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
    unit lower triangular matrix L, an n-by-n upper triangular matrix U,
    and a permutation vector piv of length m so that A(piv,:) = L*U.
    If m < n, then L is m-by-m and U is m-by-n.
    <P>
    The LU decompostion with pivoting always exists, even if the matrix is
    singular, so the constructor will never fail.  The primary use of the
    LU decomposition is in the solution of square systems of simultaneous
    linear equations.  This will fail if isNonsingular() returns false.
*/
template <typename MatrixType>
class LU
{
public :
    //typedef Real value_type;
    //typedef boost::numeric::ublas::matrix<value_type> matrix_type;
    typedef typename MatrixType::value_type value_type;
    typedef MatrixType matrix_type;
    typedef boost::numeric::ublas::vector<value_type> vector_type;
    typedef boost::numeric::ublas::vector<uint> vector_uint_type;

    /** LU Decomposition
    	@param  A   Rectangular matrix
    	@return     LU Decomposition object to access L, U and piv.
    */

    LU ( const matrix_type &A )
        :
        __LU( A.size1(), A.size2() ),
        m( A.size1() ),
        n( A.size2() ),
        pivsign( 1 ),
        piv( A.size1() )

    {
        __LU.assign(  A );

        DVLOG(2) << "LU m = "<< m << "\n";
        DVLOG(2) << "LU n = "<< n << "\n";


        // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
        for ( uint i = 0; i < m; i++ )
        {
            piv( i ) = i;
        }

        pivsign = 1;

        //value_type *LUrowi = 0;;
        vector_type LUcolj( m );

        // Outer loop.

        for ( uint j = 0; j < n; j++ )
        {

            // Make a copy of the j-th column to localize references.
            //matrix_column<matrix<value_type> > mc (m, j);
            for ( uint i = 0; i < m; i++ )
            {
                LUcolj( i ) = __LU( i,j );
            }

            // Apply previous transformations.

            for ( uint i = 0; i < m; i++ )
            {
                boost::numeric::ublas::matrix_row<matrix_type> LUrowi ( __LU, i );
                //LUrowi = __LU(i);

                // Most of the time is spent in the following dot product.

                uint kmax = std::min( i,j );
                value_type s = value_type( 0 );

                for ( uint k = 0; k < kmax; k++ )
                {
                    s += LUrowi( k )*LUcolj( k );
                }

                LUrowi( j ) = LUcolj( i ) -= s;
            }

            // Find pivot and exchange if necessary.

            uint p = j;

            for ( uint i = j+1; i < m; i++ )
            {
                if ( math::abs( LUcolj( i ) ) > math::abs( LUcolj( p ) ) )
                {
                    p = i;
                }
            }

            if ( p != j )
            {
                for ( uint k = 0; k < n; k++ )
                {
                    value_type t = __LU( p,k );
                    __LU( p,k ) = __LU( j,k );
                    __LU( j,k ) = t;
                }

                uint k = piv( p );
                piv( p ) = piv( j );
                piv( j ) = k;
                pivsign = -pivsign;
            }

            // Compute multipliers.

            if ( ( j < m ) && ( __LU( j,j ) != value_type( 0.0 ) ) )
            {
                for ( uint i = j+1; i < m; i++ )
                {
                    __LU( i,j ) /= __LU( j,j );
                }
            }
        }
    }


    /** Is the matrix nonsingular?
    	@return     1 (true)  if upper triangular factor U (and hence A)
    	is nonsingular, 0 otherwise.
    */

    uint isNonsingular ()
    {
        for ( uint j = 0; j < n; j++ )
        {
            if ( __LU( j,j ) == value_type( 0.0 ) )
                return 0;
        }

        return 1;
    }

    /**
       Return lower triangular factor
       @return     L
    */

    matrix_type getL ()
    {
        matrix_type L_( m,n );

        for ( uint i = 0; i < m; i++ )
        {
            for ( uint j = 0; j < n; j++ )
            {
                if ( i > j )
                {
                    L_( i,j ) = __LU( i,j );
                }

                else if ( i == j )
                {
                    L_( i,j ) = 1.0;
                }

                else
                {
                    L_( i,j ) = 0.0;
                }
            }
        }

        return L_;
    }

    /** Return upper triangular factor
    	@return     U portion of LU factorization.
    */

    matrix_type getU ()
    {
        matrix_type U_( n,n );

        for ( uint i = 0; i < n; i++ )
        {
            for ( uint j = 0; j < n; j++ )
            {
                if ( i <= j )
                {
                    U_( i,j ) = __LU( i,j );
                }

                else
                {
                    U_( i,j ) = 0.0;
                }
            }
        }

        return U_;
    }

    /** Return pivot permutation vector
    	@return     piv
    */
    vector_uint_type getPivot ()
    {
        return piv;
    }


    /** Compute determinant using LU factors.
    	@return     determinant of A, or 0 if A is not square.
    */
    value_type det ()
    {
        if ( m != n )
        {
            return value_type( 0 );
        }

        value_type d = value_type( pivsign );

        //DVLOG(2) << "LU::det() d= " << pivsign << "\n";
        for ( uint j = 0; j < n; j++ )
        {
            d *= __LU( j,j );
        }

        return d;
    }

    void inverse( matrix_type& __inv )
    {
        __inv = solve( ublas::identity_matrix<double>( m, m ) );
#if 0
        vector_type __t( m );
        __t = ZeroVector( __t.size() );

        for ( size_type i = 0; i < piv.size(); ++i )
        {
            __t[i] = value_type( 1 );
            ublas::row( __inv, i ) = solve( __t );
            __t[i] = value_type( 0 );
        }

#endif
    }
    /** Solve A*X = B
    	@param  B   A Matrix with as many rows as A and any number of columns.
    	@return     X so that L*U*X = B(piv,:), if B is nonconformant, returns
    	0x0 (null) array.
    */
    matrix_type solve ( const matrix_type &B )
    {

        /* Dimensions: A is mxn, X is nxk, B is mxk */

        if ( B.size1() != m )
        {
            return matrix_type( 0,0 );
        }

        if ( !isNonsingular() )
        {
            return matrix_type( 0,0 );
        }

        // Copy right hand side with pivoting
        uint nx = B.size2();


        matrix_type X ( permute_copy( B, piv, 0, nx-1 ) );

        // Solve L*Y = B(piv,:)
        for ( uint k = 0; k < n; k++ )
        {
            for ( uint i = k+1; i < n; i++ )
            {
                for ( uint j = 0; j < nx; j++ )
                {
                    X( i,j ) -= X( k,j )*__LU( i,k );
                }
            }
        }

        // Solve U*X = Y;
        for ( int k = ( int )n-1; k >= 0; k-- )
        {
            for ( uint j = 0; j < nx; j++ )
            {
                X( k,j ) /= __LU( k,k );
            }

            for ( int i = 0; i < k; i++ )
            {
                for ( uint j = 0; j < nx; j++ )
                {
                    X( i,j ) -= X( k,j )*__LU( i,k );
                }
            }
        }

        return X;
    }


    /** Solve A*x = b, where x and b are vectors of length equal
    	to the number of rows in A.

    	@param  b   a vector (Array1D> of length equal to the first dimension
    	of A.
    	@return x a vector (Array1D> so that L*U*x = b(piv), if B is nonconformant,
    	returns 0x0 (null) array.
    */

    vector_type solve ( const vector_type &b )
    {

        /* Dimensions: A is mxn, X is nxk, B is mxk */

        if ( b.size() != m )
        {
            return vector_type();
        }

        if ( !isNonsingular() )
        {
            return vector_type();
        }


        vector_type x = permute_copy( b, piv );

        // Solve L*Y = B(piv)
        for ( uint k = 0; k < n; k++ )
        {
            for ( uint i = k+1; i < n; i++ )
            {
                x( i ) -= x( k )*__LU( i,k );
            }
        }

        // Solve U*X = Y;
        for ( int k = ( int )n-1; k >= 0; k-- )
        {
            x( k ) /= __LU( k,k );

            for ( int i = 0; i < k; i++ )
            {
                x( i ) -= x( k )*__LU( i,k );
            }
        }


        return x;
    }

private:
    /* Array for internal storage of decomposition.  */
    matrix_type __LU;
    size_type m, n;
    int pivsign;
    vector_uint_type piv;


    matrix_type
    permute_copy( const matrix_type &A, const vector_uint_type &piv, uint j0, uint j1 )
    {
        uint piv_length = piv.size();

        matrix_type X( piv_length, j1-j0+1 );


        for ( uint i = 0; i < piv_length; i++ )
            for ( uint j = j0; j <= j1; j++ )
                X( i,j-j0 ) = A( piv( i ),j );

        return X;
    }

    vector_type
    permute_copy( const vector_type &A, const vector_uint_type &piv )
    {
        uint piv_length = piv.size();

        if ( piv_length != A.size() )
            return vector_type();

        vector_type x( piv_length );


        for ( uint i = 0; i < piv_length; i++ )
            x( i ) = A( piv( i ) );

        return x;
    }

}; /* class LU */

}
#endif /* __DenseLU_H */
