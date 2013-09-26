/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-03-27

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
   \file svd.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-03-27
 */
#ifndef __SVD_HPP
#define __SVD_HPP 1

#include <cmath>
#include <limits>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;
/*!
  \class SVD
  \brief Singular Value Decomposition of a rectangular matrix

  \f$ A = U * S * V^T \f$

  where matrices \f$ U\f$ and \f$ V\f$ are orthogonal and Sig is a digonal matrix.

  The singular value decomposition is performed by constructing an SVD
  object from an M*N matrix A with M>=N (that is, at least as many rows
  as columns). Note, in case M > N, matrix Sig has to be a M*N diagonal
  matrix. However, it has only N diag elements, which we store in a 1:N
  Vector sig.

  Algorithm
  Bidiagonalization with Householder reflections followed by a
  modification of a QR-algorithm. For more details, see
  G.E. Forsythe, M.A. Malcolm, C.B. Moler
  Computer methods for mathematical computations. - Prentice-Hall, 1977
  However, in the present implementation, matrices U and V are computed
  right away rather than delayed until after all Householder reflections.

  @author Christophe Prud'homme
  @see
*/
template<typename MatrixA>
class SVD
{
public:


    /** @name Typedefs
     */
    //@{

    typedef typename MatrixA::value_type value_type;
    typedef MatrixA matrix_type;
    typedef boost::numeric::ublas::vector<value_type> vector_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    SVD( matrix_type const& __A );
    ~SVD() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /** \return U of the Singular Value Decomposition */
    matrix_type const& U() const
    {
        return M_U;
    }

    /** \return V of the Singular Value Decomposition */
    matrix_type const& V() const
    {
        return M_V;
    }

    /** \return S of the Singular Value Decomposition */
    vector_type const& S() const
    {
        return M_S;
    }

    /**
       Computes \f$ C = \frac{\sigma^\mathrm{max}}{\sigma^\mathrm{min}} \f$
       which represents the condition number of the matrix that was decomposed in singular values

       \return $C$ the condition number of the matrix
    */
    value_type conditionNumber() const
    {
        return max( M_S )/min( M_S );
    }

    void conditionNumber( value_type& __max, value_type& __min ) const
    {
        __max = *std::max_element( M_S.begin(), M_S.end() );
        __min = *std::min_element( M_S.begin(), M_S.end() );
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}

protected:

private:

    SVD( SVD const & );

    /**
       \brief Left Householder Transformations

       Zero out an entire subdiagonal of the i-th column of A and compute the
       modified A[i,i] by multiplication (UP' *  A) with a matrix UP'

       (1)  UP' = E - UPi *  UPi' / beta

       where a column-vector UPi is as follows

       (2)  UPi = [ (i-1) zeros, A[i,i] + Norm, vector A[i+1:M,i] ]

       where beta = UPi'   A[,i] and Norm is the norm of a vector A[i:M,i]
       (sub-diag part of the i-th column of A). Note we assign the Norm the
       same sign as that of A[i,i].

       By construction, (1) does not affect the first (i-1) rows of A. Since
       A[*,1:i-1] is bidiagonal (the result of the i-1 previous steps of
       the bidiag algorithm), transform (1) doesn't affect these i-1 columns
       either as one can easily verify.

       The i-th column of A is transformed as

       (3)  UP'   A[*,i] = A[*,i] - UPi

       (since UPi'*A[*,i]/beta = 1 by construction of UPi and beta)

       This means effectively zeroing out A[i+1:M,i] (the entire subdiagonal
       of the i-th column of A) and replacing A[i,i] with the -Norm. Thus
       modified A[i,i] is returned by the present function.

       The other (i+1:N) columns of A are transformed as

       (4)  UP'   A[,j] = A[,j] - UPi * ( UPi' * A[,j] / beta )

       Note, due to (2), only elements of rows i+1:M actually  participate
       in above transforms; the upper i-1 rows of A are not affected.
       As was mentioned earlier,

       (5)  beta = UPi' * A[,i] = (A[i,i] + Norm)*A[i,i] + A[i+1:M,i]*A[i+1:M,i]
       = ||A[i:M,i]||^2 + Norm*A[i,i] = Norm^2 + Norm*A[i,i]

       (note the sign of the Norm is the same as A[i,i])
       For extra precision, vector UPi (and so is Norm and beta) are scaled,
       which would not affect (4) as easy to see.

       To satisfy the definition

       (6)  .SIG = U' * A * V

       the result of consecutive transformations (1) over matrix A is accumulated
       in matrix U' (which is initialized to be a unit matrix). At each step,
       U' is left-multiplied by UP' = UP (UP is symmetric by construction,
       see (1)). That is, U is right-multiplied by UP, that is, rows of U are
       transformed similarly to columns of A, see eq. (4). We also keep in mind
       that multiplication by UP at the i-th step does not affect the first i-1
       columns of U.

       Note that the vector UPi doesn't have to be allocated explicitly: its
       first i-1 components are zeros (which we can always imply in computations),
       and the rest of the components (but the UPi[i]) are the same as those
       of A[i:M,i], the subdiagonal of A[,i]. This column, A[,i] is affected only
       trivially as explained above, that is, we don't need to carry this
       transformation explicitly (only A[i,i] is going to be non-trivially
       affected, that is, replaced by -Norm, but we will use sig[i] to store
       the result).

    */
    value_type leftHouseholder( MatrixA& A, const int i );


    /**
       \brief Right Householder Transformations

       Zero out i+2:N elements of a row A[i,] of matrix A by right
       multiplication (A   VP) with a matrix VP
       (1)  VP = E - VPi   VPi' / beta

       where a vector-column .VPi is as follows

       (2)  VPi = [ i zeros, A[i,i+1] + Norm, vector A[i,i+2:N] ]

       where beta = A[i,]   VPi and Norm is the norm of a vector A[i,i+1:N]
       (right-diag part of the i-th row of A). Note we assign the Norm the
       same sign as that of A[i,i+1].

       By construction, (1) does not affect the first i columns of A. Since
       A[1:i-1,] is bidiagonal (the result of the previous steps of
       the bidiag algorithm), transform (1) doesn't affect these i-1 rows
       either as one can easily verify.

       The i-th row of A is transformed as

       (3)  A[i,*]   VP = A[i,*] - VPi'

       (since A[i,*]*VPi/beta = 1 by construction of VPi and beta)
       This means effectively zeroing out A[i,i+2:N] (the entire right super-
       diagonal of the i-th row of A, but ONE superdiag element) and replacing
       A[i,i+1] with - Norm. Thus modified A[i,i+1] is returned as the result of
       the present function.

       The other (i+1:M) rows of A are transformed as

       (4)  A[j,]   VP = A[j,] - VPi'   ( A[j,]   VPi / beta )

       Note, due to (2), only elements of columns i+1:N actually  participate
       in above transforms; the left i columns of A are not affected.

       As was mentioned earlier,

       (5)  beta = A[i,]   VPi = (A[i,i+1] + Norm)*A[i,i+1]+ A[i,i+2:N]*A[i,i+2:N]
       = ||A[i,i+1:N]||^2 + Norm*A[i,i+1] = Norm^2 + Norm*A[i,i+1]

       (note the sign of the Norm is the same as A[i,i+1])
       For extra precision, vector VPi (and so is Norm and beta) are scaled,
       which would not affect (4) as easy to see.

       The result of consecutive transformations (1) over matrix A is accumulated
       in matrix V (which is initialized to be a unit matrix). At each step,
       V is right-multiplied by VP. That is, rows of V are transformed similarly
       to rows of A, see eq. (4). We also keep in mind that multiplication by
       VP at the i-th step does not affect the first i rows of V.
       Note that vector VPi doesn't have to be allocated explicitly: its
       first i components are zeros (which we can always imply in computations),
       and the rest of the components (but the VPi[i+1]) are the same as those
       of A[i,i+1:N], the superdiagonal of A[i,]. This row, A[i,] is affected
       only trivially as explained above, that is, we don't need to carry this
       transformation explicitly (only A[i,i+1] is going to be non-trivially
       affected, that is, replaced by -Norm, but we will use super_diag[i+1] to
       store the result).

    */
    value_type rightHouseholder( MatrixA& A, const int i );

    /**
       \brief Bidiagonalization
       This nethod turns matrix A into a bidiagonal one. Its N diagonal elements
       are stored in a vector sig, while N-1 superdiagonal elements are stored
       in a vector super_diag(2:N) (with super_diag(1) being always 0).
       Matrices U and V store the record of orthogonal Householder
       reflections that were used to convert A to this form. The method
       returns the norm of the resulting bidiagonal matrix, that is, the
       maximal column sum.
    */
    value_type bidiagonalize( vector_type& super_diag, const matrix_type& _A );

    /*
      \brief QR-diagonalization of a bidiagonal matrix

      After bidiagonalization we get a bidiagonal matrix J:
      (1)  J = U' * A * V
      The present method turns J into a matrix JJ by applying a set of
      orthogonal transforms
      (2)  JJ = S' * J * T
      Orthogonal matrices S and T are chosen so that JJ were also a
      bidiagonal matrix, but with superdiag elements smaller than those of J.
      We repeat (2) until non-diag elements of JJ become smaller than EPS
      and can be disregarded.
      Matrices S and T are constructed as
      (3)  S = S1 * S2 * S3 ... Sn, and similarly T
      where Sk and Tk are matrices of simple rotations
      (4)  Sk[i,j] = i==j ? 1 : 0 for all i>k or i<k-1
      Sk[k-1,k-1] = cos(Phk),  Sk[k-1,k] = -sin(Phk),
      SK[k,k-1] = sin(Phk),    Sk[k,k] = cos(Phk), k=2..N
      Matrix Tk is constructed similarly, only with angle Thk rather than Phk.

      Thus left multiplication of J by SK' can be spelled out as
      (5)  (Sk' * J)[i,j] = J[i,j] when i>k or i<k-1,
      [k-1,j] = cos(Phk)*J[k-1,j] + sin(Phk)*J[k,j]
      [k,j] =  -sin(Phk)*J[k-1,j] + cos(Phk)*J[k,j]
      That is, k-1 and k rows of J are replaced by their linear combinations;
      the rest of J is unaffected. Right multiplication of J by Tk similarly
      changes only k-1 and k columns of J.
      Matrix T2 is chosen the way that T2'J'JT2 were a QR-transform with a
      shift. Note that multiplying J by T2 gives rise to a J[2,1] element of
      the product J (which is below the main diagonal). Angle Ph2 is then
      chosen so that multiplication by S2' (which combines 1 and 2 rows of J)
      gets rid of that elemnent. But this will create a [1,3] non-zero element.
      T3 is made to make it disappear, but this leads to [3,2], etc.
      In the end, Sn removes a [N,N-1] element of J and matrix S'JT becomes
      bidiagonal again. However, because of a special choice
      of T2 (QR-algorithm), its non-diag elements are smaller than those of J.

      All this process in more detail is described in
      J.H. Wilkinson, C. Reinsch. Linear algebra - Springer-Verlag, 1971

      If during transforms (1), JJ[N-1,N] turns 0, then JJ[N,N] is a singular
      number (possibly with a wrong (that is, negative) sign). This is a
      consequence of Frantsis' Theorem, see the book above. In that case, we can
      eliminate the N-th row and column of JJ and carry out further transforms
      with a smaller matrix. If any other superdiag element of JJ turns 0,
      then JJ effectively falls into two independent matrices. We will process
      them independently (the bottom one first).

      Since matrix J is a bidiagonal, it can be stored efficiently. As a matter
      of fact, its N diagonal elements are in array Sig, and superdiag elements
      are stored in array super_diag.

      Carry out U * S with a rotation matrix
      S (which combines i-th and j-th columns
      of U, i>j)

    */
    void rotate( matrix_type& U, const int i, const int j, const value_type cos_ph, const value_type sin_ph );

    /*
      A diagonal element J[l-1,l-1] turns out 0 at the k-th step of the
      algorithm. That means that one of the original matrix' singular numbers
      shall be zero. In that case, we multiply J by specially selected
      matrices S' of simple rotations to eliminate a superdiag element J[l-1,l].
      After that, matrix J falls into two pieces, which can be dealt with
      in a regular way (the bottom piece first).

      These special S transformations are accumulated into matrix U: since J
      is left-multiplied by S', U would be right-multiplied by S. Transform
      formulas for doing these rotations are similar to (5) above. See the
      book cited above for more details.
    */
    void ripThrough( vector_type& super_diag, const int k, const int l, const double eps );

    /**
       We're about to proceed doing QR-transforms
       on a (bidiag) matrix J[1:k,1:k]. It may happen
       though that the matrix splits (or can be
       split) into two independent pieces. This function
       checks for splitting and returns the lowerbound
       index l of the bottom piece, J[l:k,l:k]
    */
    int getWorkSubmatrix( vector_type& super_diag, const int k, const double eps );

    /** main algorithm */
    void diagonalize( vector_type& super_diag, const double eps );

private:

    uint M_rows;
    uint M_cols;

    matrix_type M_U;
    matrix_type M_V;
    vector_type M_S;

};

template<typename MatrixA>
class SOrth
{
public:

    /** @name Typedefs
     */
    //@{

    typedef typename MatrixA::value_type value_type;
    typedef MatrixA matrix_type;
    typedef boost::numeric::ublas::vector<value_type> vector_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    SOrth( matrix_type const& __A );
    ~SOrth() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /** \return U of the Singular Value Decomposition */
    matrix_type const& Q() const
    {
        return M_Q;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{


    //@}
private:

    uint M_rows;
    uint M_cols;
    matrix_type M_Q;
};

/**
   \todo unfinished implementation: need to extract
*/
template<typename MatrixA>
SOrth<MatrixA>::SOrth( MatrixA const& A )
    :
    M_rows( A.size1() ),
    M_cols( A.size2() ),
    M_Q()
{
    using namespace boost::numeric::ublas;

    SVD<matrix_type> __svd( A );
    vector_type __S;

    if ( M_rows > 1 )
        __S = __svd.S();

    else if ( M_rows == 1 )
        __S( 0 ) = __svd.S()( 0 );

    const value_type eps = std::numeric_limits<value_type>::epsilon() * *std::max_element( __S.begin(), __S.end() )* std::max( M_rows, M_cols );

    /*[m,n] = size(A);
      if m > 1, s = diag(S);
      elseif m == 1, s = S(1);
      else s = 0;
      end
      tol = max(m,n) * max(s) * eps;
      r = sum(s > tol);
      Q = U(:,1:r);*/

}
/*
 * SVD
 */

template<typename MatrixA>
SVD<MatrixA>::SVD( MatrixA const& A )
    :
    M_rows( A.size1() ),
    M_cols( A.size2() ),
    M_U( M_rows, M_rows ),
    M_V( M_cols, M_cols ),
    M_S( M_cols )
{
    //FEELPP_ASSERT( M_rows >= M_cols ).error("The Matrix A should have at least as many rows as it has columns");

    MatrixA B( A );

    bool do_transpose = M_rows < M_cols;

    if ( do_transpose )
    {
        // transpose
        M_rows = A.size2();
        M_cols = A.size1();
        M_U.resize( M_rows, M_rows );
        M_V.resize( M_cols, M_cols );
        M_S.resize( M_cols );
        B = ublas::trans( A );
    }

    M_U = boost::numeric::ublas::identity_matrix<value_type>( M_rows );
    M_V = boost::numeric::ublas::identity_matrix<value_type>( M_cols );

    vector_type super_diag( M_cols );

    const value_type bidiag_norm = bidiagonalize( super_diag,B );

    // define Significance threshold
    const value_type eps = std::numeric_limits<value_type>::epsilon() * bidiag_norm;

    diagonalize( super_diag,eps );

    if ( do_transpose )
    {
        matrix_type Temp( M_V );
        M_V = ublas::trans( M_U );
        M_U = ublas::trans( Temp );
    }
}

/*
  Bidiagonalization
*/


template<typename MatrixA>
typename SVD<MatrixA>::value_type
SVD<MatrixA>::leftHouseholder( matrix_type& A, const int i )
{
    value_type scale = 0;                 // Compute the scaling factor

    for ( int k=i; k< M_rows; k++ )
        scale += math::abs( A( k,i ) );

    if ( scale == 0 )                    // If A[,i] is a null vector, no
        return 0;                          // transform is required

    value_type Norm_sqr = 0;                 // Scale UPi (that is, A[,i])

    for ( int k=i; k< M_rows; k++ )                // and compute its norm, Norm^2
    {
        A( k,i ) /= scale;
        Norm_sqr +=  A( k,i )*A( k,i );
    }

    value_type new_Aii = math::sqrt( Norm_sqr );   // new_Aii = -Norm, Norm has the

    if ( A( i,i ) > 0 ) new_Aii = -new_Aii; // same sign as Aii (that is, UPi[i])

    const value_type beta = - A( i,i )*new_Aii + Norm_sqr;
    A( i,i ) -= new_Aii;                 // UPi[i] = A[i,i] - (-Norm)

    for ( int j=i+1; j<M_cols; j++ )      // Transform i+1:N columns of A
    {
        value_type factor = 0;

        for ( int k=i; k< M_rows; k++ )
            factor += A( k,i ) * A( k,j );

        factor /= beta;

        for ( int k=i; k< M_rows; k++ )
            A( k,j ) -= A( k,i ) * factor;
    }

    for ( size_type j=0; j< M_rows; j++ )                // Accumulate the transform in U
    {
        value_type factor = 0;

        for ( size_type k=i; k< M_rows; k++ )
            factor += A( k,i ) * M_U( j,k );  // Compute  U[j,] * UPi

        factor /= beta;

        for ( size_type k=i; k< M_rows; k++ )
            M_U( j,k ) -= A( k,i ) * factor;
    }

    return new_Aii * scale;              // Scale new Aii back (our new Sig[i])
}

template<typename MatrixA>
typename SVD<MatrixA>::value_type
SVD<MatrixA>::rightHouseholder( matrix_type& A, const int i )
{
    value_type scale = 0;                 // Compute the scaling factor

    for ( uint16_type k=i+1; k<M_cols; k++ )
        scale += math::abs( A( i,k ) );

    if ( scale == 0 )                    // If A[i,] is a null vector, no
        return 0;                          // transform is required

    value_type Norm_sqr = 0;                 // Scale VPi (that is, A[i,])

    for ( uint16_type k=i+1; k<M_cols; k++ )      // and compute its norm, Norm^2
    {
        A( i,k ) /= scale;
        Norm_sqr += A( i,k )*A( i,k );
    }

    value_type new_Aii1 = math::sqrt( Norm_sqr );  // new_Aii1 = -Norm, Norm has the

    if ( A( i,i+1 ) > 0 )                // same sign as
        new_Aii1 = -new_Aii1;              // Aii1 (that is, VPi[i+1])

    const value_type beta = - A( i,i+1 )*new_Aii1 + Norm_sqr;
    A( i,i+1 ) -= new_Aii1;              // VPi[i+1] = A[i,i+1] - (-Norm)

    for ( uint16_type j=i+1; j<M_rows; j++ )      // Transform i+1:M rows of A
    {
        value_type factor = 0;

        for ( uint16_type k=i+1; k<M_cols; k++ )
            factor += A( i,k ) * A( j,k );   // Compute A[j,] * VPi

        factor /= beta;

        for ( uint16_type k=i+1; k<M_cols; k++ )
            A( j,k ) -= A( i,k ) * factor;
    }

    for ( uint16_type j=0; j<M_cols; j++ )                // Accumulate the transform in V
    {
        value_type factor = 0;

        for ( uint16_type k=i+1; k<M_cols; k++ )
            factor += A( i,k ) * M_V( j,k );  // Compute  V[j,] * VPi

        factor /= beta;

        for ( uint16_type k=i+1; k<M_cols; k++ )
            M_V( j,k ) -= A( i,k ) * factor;
    }

    return new_Aii1 * scale;             // Scale new Aii1 back
}

template<typename MatrixA>
typename SVD<MatrixA>::value_type
SVD<MatrixA>::bidiagonalize( vector_type& super_diag, const matrix_type& _A )
{
    value_type __norm_acc = 0;

    // No superdiag elem above A(0,0)
    super_diag( 0 ) = 0;

    // A being transformed
    matrix_type A = _A;

    for ( uint16_type __i = 0; __i < M_cols; ++__i )
    {
        M_S( __i ) = leftHouseholder( A,__i );
        //std::cout << "S=" << M_S << "\n";
        const value_type& diagi = M_S( __i );

        if ( __i < M_cols-1 )
            super_diag( __i+1 ) = rightHouseholder( A,__i );

        __norm_acc = std::max( __norm_acc, math::abs( diagi ) + math::abs( super_diag( __i ) ) );
    }

    return __norm_acc;
}


template<typename MatrixA>
void
SVD<MatrixA>::rotate( matrix_type& U, const int i, const int j,
                      const value_type cos_ph,  const value_type sin_ph )
{
    using namespace boost::numeric::ublas;

    matrix_column<matrix_type> Ui ( U, i );
    matrix_column<matrix_type> Uj ( U, j );

    for ( unsigned __k = 0; __k < Ui.size(); ++ __k )
    {
        const value_type Ujk_was = Uj( __k );

        Uj( __k ) = cos_ph  * Ujk_was + sin_ph * Ui( __k );
        Ui( __k ) = -sin_ph * Ujk_was + cos_ph * Ui( __k );
    }
}


template<typename MatrixA>
void
SVD<MatrixA>::ripThrough( vector_type& super_diag, const int k, const int l, const double eps )
{
    // Accumulate cos,sin of Ph
    value_type cos_ph = 0, sin_ph = 1;

    // The first step of the loop below
    // (when i==l) would eliminate J[l-1,l],
    // which is stored in super_diag(l)
    // However, it gives rise to J[l-1,l+1]
    // and J[l,l+2]
    // The following steps eliminate these
    // until they fall below
    // significance
    for (  int i=l; i<=k; i++ )
    {
        const value_type f = sin_ph * super_diag( i );

        super_diag( i ) *= cos_ph;

        // Current J[l-1,l] may become unsignificant
        if ( math::abs( f ) <= eps )
            break;

        // unnormalized sin/cos
        cos_ph = M_S( i );
        sin_ph = -f;

        // sqrt(sin^2+cos^2)
        const value_type norm = ( M_S( i ) = hypot( cos_ph,sin_ph ) );

        // Normalize sin/cos
        cos_ph /= norm;
        sin_ph /= norm;

        rotate( M_U, i, l-1, cos_ph, sin_ph );
    }
}
template<typename MatrixA>
int
SVD<MatrixA>::getWorkSubmatrix( vector_type& super_diag, const int k, const double eps )
{
    for ( int l=k; l > 0; --l )
    {
        if ( math::abs( super_diag( l ) ) <= eps )
        {
            // The breaking point: zero J[l-1,l]
            return l;
        }

        else if ( math::abs( M_S( l-1 ) ) <= eps )
        {
            // Diagonal J[l,l] turns out 0
            // meaning J[l-1,l] _can_ be made
            // zero after some rotations
            ripThrough( super_diag,k,l,eps );
            return l;
        }
    }

    // Deal with J[1:k,1:k] as a whole
    return 0;
}

template<typename MatrixA>
void
SVD<MatrixA>::diagonalize( vector_type& super_diag, const double eps )
{
    for ( int k=M_cols-1; k >= 0; --k )	// QR-iterate upon J[l:k,l:k]
    {
        //std::cerr << "k = " << k << "\n";
        int l;

        // until superdiag J[k-1,k] becomes 0
        while ( math::abs( super_diag( k ) ) > eps )
        {
            l=getWorkSubmatrix( super_diag,k,eps );
            //std::cerr << "l = " << l << "\n";
            value_type shift;			// Compute a QR-shift from a bottom
            {
                // corner minor of J[l:k,l:k] order 2
                value_type Jk2k1 = super_diag( k-1 ),	// J[k-2,k-1]
                           Jk1k  = super_diag( k ),
                           Jk1k1 = M_S( k-1 ),		// J[k-1,k-1]
                           Jkk   = M_S( k ),
                           Jll   = M_S( l );		// J[l,l]
                shift = ( Jk1k1-Jkk )*( Jk1k1+Jkk ) + ( Jk2k1-Jk1k )*( Jk2k1+Jk1k );
                shift /= 2*Jk1k*Jk1k1;
                shift += ( shift < 0 ? -1 : 1 ) * math::sqrt( shift*shift+1 );
                shift = ( ( Jll-Jkk )*( Jll+Jkk ) + Jk1k*( Jk1k1/shift-Jk1k ) )/Jll;
            }
            // Carry on multiplications by T2, S2, T3...
            value_type cos_th = 1;
            value_type sin_th = 1;
            value_type Ji1i1 = M_S( l );	// J[i-1,i-1] at i=l+1...k

            for ( int i=l+1; i<=k; i++ )
            {
                value_type Ji1i = super_diag( i ); // J[i-1,i]
                value_type Jii = M_S( i ); //  J[i,i]

                sin_th *= Ji1i;
                Ji1i *= cos_th;
                cos_th = shift;

                value_type norm_f = ( super_diag( i-1 ) = hypot( cos_th,sin_th ) );
                cos_th /= norm_f, sin_th /= norm_f;

                // Rotate J[i-1:i,i-1:i] by Ti
                shift = cos_th*Ji1i1 + sin_th*Ji1i;	// new J[i-1,i-1]
                Ji1i = -sin_th*Ji1i1 + cos_th*Ji1i;	// J[i-1,i] after rotation
                const value_type Jii1 = Jii*sin_th;		// Emerged J[i,i-1]
                Jii *= cos_th;				// new J[i,i]

                // Accumulate T rotations in V
                rotate( M_V, i, i-1, cos_th, sin_th );

                value_type cos_ph = shift, sin_ph = Jii1;// Make Si to get rid of J[i,i-1]
                M_S( i-1 ) = ( norm_f = hypot( cos_ph,sin_ph ) );	// New J[i-1,i-1]

                if ( norm_f == 0 )		// If norm =0, rotation angle
                {
                    cos_ph = cos_th;
                    sin_ph = sin_th; // can be anything now
                }

                else
                {
                    cos_ph /= norm_f;
                    sin_ph /= norm_f;
                }

                // Rotate J[i-1:i,i-1:i] by Si
                shift = cos_ph * Ji1i + sin_ph*Jii;	// New J[i-1,i]
                Ji1i1 = -sin_ph*Ji1i + cos_ph*Jii;	// New Jii, would carry over as J[i-1,i-1] for next i

                // Accumulate S rotations in U
                rotate( M_U, i, i-1, cos_ph, sin_ph );

                // Jii1 disappears, sin_th would
                cos_th = cos_ph, sin_th = sin_ph; // carry over a (scaled) J[i-1,i+1]
                // to eliminate on the next i, cos_th
                // would carry over a scaled J[i,i+1]
            }

            super_diag( l ) = 0;		// Supposed to be eliminated by now
            super_diag( k ) = shift;
            M_S( k ) = Ji1i1;
        }		// --- end-of-QR-iterations

        if ( M_S( k ) < 0 )		// Correct the sign of the sing number
        {
            M_S( k ) = -M_S( k );

            using namespace boost::numeric::ublas;
            //matrix_column<matrix<value_type> > Vk ( M_V, k);
            matrix_column<matrix_type> Vk ( M_V, k );

            for ( unsigned l = 0; l < Vk.size(); ++ l )
            {
                Vk( l ) = -Vk( l );
            }
        }
    }
}


}
#endif  /* __SVD_HPP  */
