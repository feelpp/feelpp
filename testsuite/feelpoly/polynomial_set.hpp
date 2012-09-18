/**
   \file Polynomial_set.hpp
   \author Abdoulaye Samake <samakeablo@gmail.com>
   \date 2010-05-28
 */

/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-07-02

  Copyright (C) 2010 Université Joseph Fourier Grenoble 1

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
   \file Polynomial_set.hpp
   \author Abdoulaye Samake <samakeablo@gmail.com>
   \date 2010-07-02
 */


# ifndef DEF_Polynomial_set_H
# define DEF_Polynomial_set_H

#include <vector>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>



namespace Feel
{


/**
 * \class Polynomial_set
 * \brief polynomial_set class
 *
 * The polynomial is expressed in the basis from \p Poly. The
 * coefficients of the polynomial in this basis are represented by a
 * matrix whose lines are the polymomial components coefficients (1 if
 * \code is_scalar == true \endcode, \p nDim if \code is_vectorial ==
 * true\endcode and columns are the basis
 *
 * Evaluating the polynomial at a set of points(or just one point) is
 * then simply a matrix-matrix product.
 *
 * \ingroup Polynomial_set
 * @author Abdoulaye Samake < samakeablo@gmail.com>
 * @see
 */
template< uint nDim ,uint nOrder , uint dim ,uint number >
class Polynomialset
{

public:

    /** @name Typedefs
     */
    //@{

    typedef  boost::numeric::ublas::matrix<double> matrix_type ;
    typedef  boost::numeric::ublas::vector<double> vector_type ;
    typedef  boost::numeric::ublas::identity_matrix<double> matrix_id ;
    typedef  boost::numeric::ublas::zero_vector<double> vector_zero ;
    typedef  boost::numeric::ublas::vector<double> Points_type ;
    typedef  boost::numeric::ublas::vector<Points_type> vectors_type ;

    //@}

    //@}

    /** @name Constructors, destructor
     */
    //@

    /**
     * default constructor
     */

    Polynomialset() ;


    Polynomialset ( matrix_type const&  m , vector_type const&  b ) ;
    Polynomialset ( matrix_type const&  m ) ;
    Polynomialset ( const Polynomialset<nDim,nOrder ,dim ,number> & P ) ;
    ~Polynomialset() ;
    matrix_type evaluates_Points( vectors_type const&  P );

private :

    matrix_type M ;

    vector_type B ;

}; // Polynomial_set

// Implementation

template< uint nDim ,uint nOrder , uint dim ,uint number >
Polynomialset <nDim , nOrder ,dim ,number>::Polynomialset() :M( matrix_id( nOrder,dim ) ),B( vector_zero( nDim ) )
{
};

template< uint nDim ,uint nOrder ,uint dim ,uint number >
Polynomialset <nDim , nOrder , dim ,number>::Polynomialset ( matrix_type const&  m , vector_type const&  b )
    :
    M( m ),
    B(  b )

{
}




template< uint nDim , uint nOrder ,uint dim ,uint number >
Polynomialset <nDim , nOrder , dim ,number>::Polynomialset ( matrix_type const&  m )
    :
    M( m ),
    B( ublas::scalar_vector<double>( nDim, 1. ) )
{
}




template< uint nDim , uint nOrder ,uint dim ,uint number >
Polynomialset<nDim ,nOrder ,dim ,number >::Polynomialset ( const Polynomialset< nDim, nOrder ,dim ,number> & P )
    :
    M( P.M ),
    B( P.B )

{
}



template< uint nDim , uint nOrder ,uint dim ,uint number >
Polynomialset<nDim ,nOrder ,dim ,number>::~Polynomialset()
{
}

template< uint nDim , uint nOrder ,uint dim ,uint number >
typename Polynomialset<nDim ,nOrder ,dim ,number>::matrix_type
Polynomialset<nDim ,nOrder ,dim ,number>::evaluates_Points( vectors_type const&  P )
{
    matrix_type matrice( nDim , number ) ;
    matrix_type matrice_valeur( nDim , number ) ;


    if ( P( 0 ).size() == 2 )
    {
        for ( uint j = 0 ; j< matrice.size2() ; ++j )
        {
            matrice( 0 ,j ) =  1. ;
            matrice( 1 ,j ) = B( 1 )*P( j )( 0 ) ;
            matrice( 2 ,j ) = B( 2 )*P( j )( 1 ) ;
            matrice( 3 ,j ) = B( 3 )*P( j )( 0 )*P( j )( 1 ) ;
            matrice( 4 ,j ) = B( 4 )*pow( P( j )( 0 ),2 ) ;
            matrice( 5 ,j ) = B( 5 )*pow( P( j )( 1 ),2 ) ;
            matrice( 6 ,j ) = B( 6 )*pow( P( j )( 0 ),2 )*P( j )( 1 ) ;
            matrice( 7 ,j ) = B( 7 )*pow( P( j )( 1 ),2 )*P( j )( 0 ) ;
            matrice( 8 ,j ) = B( 8 )*pow( P( j )( 0 ),3 ) ;
            matrice( 9 ,j ) = B( 9 )*pow( P( j )( 1 ),3 ) ;

        }
    }

    else
    {
        for ( uint j = 0 ; j< matrice.size2() ; ++j )
        {
            matrice( 0 ,j ) =  1. ;
            matrice( 1 ,j ) = B( 1 )*P( j )( 0 ) ;
            matrice( 2 ,j ) = B( 2 )*P( j )( 1 ) ;
            matrice( 3 ,j ) = B( 3 )*P( j )( 2 ) ;
            matrice( 4 ,j ) = B( 4 )*P( j )( 0 )*P( j )( 1 ) ;
            matrice( 5 ,j ) = B( 5 )*P( j )( 0 )*P( j )( 2 ) ;
            matrice( 6 ,j ) = B( 6 )*P( j )( 1 )*P( j )( 2 ) ;
            matrice( 7 ,j ) = B( 7 )*pow( P( j )( 0 ),2 ) ;
            matrice( 8 ,j ) = B( 8 )*pow( P( j )( 1 ),2 ) ;
            matrice( 9 ,j ) = B( 9 )*pow( P( j )( 2 ),2 ) ;
            matrice( 10 ,j ) = B( 10 )*pow( P( j )( 0 ),2 )*P( j )( 1 ) ;
            matrice( 11 ,j ) = B( 11 )*pow( P( j )( 1 ),2 )*P( j )( 0 ) ;
            matrice( 12 ,j ) = B( 12 )*pow( P( j )( 0 ),2 )*P( j )( 2 ) ;
            matrice( 13 ,j ) = B( 13 )*pow( P( j )( 2 ),2 )*P( j )( 0 ) ;
            matrice( 14 ,j ) = B( 14 )*pow( P( j )( 1 ),2 )*P( j )( 2 ) ;
            matrice( 15 ,j ) = B( 15 )*pow( P( j )( 2 ),2 )*P( j )( 1 ) ;
            matrice( 16 ,j ) = B( 16 )*P( j )( 0 )*P( j )( 1 )*P( j )( 2 ) ;
            matrice( 17 ,j ) = B( 17 )*pow( P( j )( 0 ),3 ) ;
            matrice( 18 ,j ) = B( 18 )*pow( P( j )( 1 ),3 ) ;
            matrice( 19 ,j ) = B( 19 )*pow( P( j )( 2 ),3 ) ;

        }
    }

#if 0

    for ( uint i = 0 ; i < matrice_valeur.size1() ; ++i )
    {
        for ( uint j = 0 ; j< matrice_valeur.size2() ; ++j )
        {
            matrice_valeur ( i , j ) = matrice ( i , j ) ;
        }
    }

#endif
    return   boost::numeric::ublas::prod( M , matrice )  ;
}



}  //namespace Feel

#endif









