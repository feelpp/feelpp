/* -*- mode: c++ -*-

  This file is part of the Life library

 Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
                 Abdoulaye Samake <samakeablo@gmail.com>
       Date: 2010-05-28

 Copyright (C) 2010 Université Joseph Fourier Grenoble 1

*/
/**
   \file Polynomial_set.hpp
   \author Abdoulaye Samake <samakeablo@gmail.com>
   \date 2010-05-28
 */


# ifndef DEF_Polynomial_set_H
# define DEF_Polynomial_set_H

#include <vector>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>


namespace Life
{
 typedef  boost::numeric::ublas::matrix<double> matrix_type ;
 typedef  boost::numeric::ublas::vector<double> vector_type ;
 typedef  boost::numeric::ublas::identity_matrix<double> matrix_id ;
 typedef  boost::numeric::ublas::zero_vector<double> vector_zero ;
 typedef  boost::numeric::ublas::vector<double> Points_type ;
 typedef  boost::numeric::ublas::vector<Points_type> vectors_type ;


template< uint nDim ,uint nOrder , uint dim ,uint number >
class Polynomialset
{

public:



    Polynomialset() ;
    Polynomialset ( matrix_type  m , vector_type  b) ;
    Polynomialset ( matrix_type  m ) ;
    Polynomialset ( const Polynomialset<nDim,nOrder ,dim ,number> & P) ;
   ~Polynomialset() ;
    matrix_type evaluate_Points( vectors_type P ) ;

private :

    matrix_type M ;

    vector_type B ;

 }; // Polynomial_set

// Implementation

template< uint nDim ,uint nOrder , uint dim ,uint number >
Polynomialset <nDim , nOrder ,dim ,number>::Polynomialset() :M( matrix_id(nOrder,dim)),B(vector_zero(nDim))
{
};

template< uint nDim ,uint nOrder ,uint dim ,uint number >
Polynomialset <nDim , nOrder , dim ,number>::Polynomialset ( matrix_type  m , vector_type  b)
{
    M =  matrix_type(nOrder,dim) ;

    B =  vector_type(nDim) ;
    for(uint i = 0 ;i < M.size1() ; ++i)
    {
       for(uint j = 0 ;j< M.size2() ; ++j)
       {
           M(i , j) = m(i , j) ;
       }
    }
  for(uint k = 0 ;k< B.size() ; ++k)
    {
      B( k ) = b( k ) ;
    }


}




template< uint nDim , uint nOrder ,uint dim ,uint number >
Polynomialset <nDim , nOrder , dim ,number>::Polynomialset ( matrix_type  m )
{
    M =  matrix_type(nOrder,dim) ;
    B =  vector_type(nDim) ;
    for(uint i = 0 ;i < M.size1() ; ++i)
    {
        for(uint j = 0 ;j< M.size2() ; ++j)
        {
            M(i , j) = m(i , j) ;
        }
    }
    for(uint k = 0 ;k< B.size() ; ++k)
    {
        B( k ) = 1. ;
    }

}




template< uint nDim , uint nOrder ,uint dim ,uint number >
Polynomialset<nDim ,nOrder ,dim ,number >::Polynomialset ( const Polynomialset< nDim, nOrder ,dim ,number> & P)
{
    M =  matrix_type(nOrder,dim) ;

    B =  vector_type(nDim) ;
    for(uint i = 0 ;i < M.size1() ; ++i)
    {
        for(uint j = 0 ;j< M.size2() ; ++j)
        {
            M(i , j) = P.M(i , j) ;
        }
    }


    for(uint k = 0 ;k< B.size() ; ++k)
    {
        B( k ) = P.B( k ) ;
    }

}



template< uint nDim , uint nOrder ,uint dim ,uint number >
Polynomialset<nDim ,nOrder ,dim ,number>::~Polynomialset()
{
}

template< uint nDim , uint nOrder ,uint dim ,uint number >
matrix_type Polynomialset<nDim ,nOrder ,dim ,number>::evaluate_Points( vectors_type  P)
{
  matrix_type matrice(nDim , number) ;
  matrix_type matrice_valeur(dim , number) ;



    for(uint j = 0 ;j< matrice.size2() ;++j)
    {
       matrice(0 ,j) =  1. ;
       matrice(1 ,j) = B(1)*P(j)(0) ;
       matrice(2 ,j) = B(2)*P(j)(1) ;
       matrice(3 ,j) = B(3)*P(j)(0)*P(j)(1) ;
       matrice(4 ,j) = B(4)*pow(P(j)(0),2) ;
       matrice(5 ,j) = B(5)*pow(P(j)(1),2) ;
       matrice(6 ,j) = B(6)*pow(P(j)(0),2)*P(j)(1) ;
       matrice(7 ,j) = B(7)*pow(P(j)(1),2)*P(j)(0) ;
       matrice(8 ,j) = B(8)*pow(P(j)(0),3) ;
       matrice(9 ,j) = B(9)*pow(P(j)(1),2) ;
    }

    for(uint i = 0 ;i < matrice_valeur.size1() ;++i)
  {
    for(uint j = 0 ;j< matrice_valeur.size2() ;++j)
    {
        matrice_valeur ( i , j ) = matrice ( i , j ) ;
    }
  }



 return   boost::numeric::ublas::prod( M , matrice_valeur )  ;




}




}  //namespace Life

#endif









