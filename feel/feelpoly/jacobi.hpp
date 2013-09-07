/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-26

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2012 Universite Joseph Fourier Grenoble 1

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
   \file jacobi.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-26
 */
#ifndef __Jacobi_H
#define __Jacobi_H 1

#include <cmath>

#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#if 0
#include <Eigen/Core>
#endif

#include <feel/feelcore/traits.hpp>

namespace Feel
{
namespace ublas=boost::numeric::ublas;
/**
 * \class Jacobi
 * \brief 1D Jacobi polynomial
 *
 * (excerpt from Karniadakis/Sherwin Appendix A)
 * Jacobi polynomials \f$ P^{\alpha,beta}_n(x)\f$ are a family of
 * polynomial to the singular Sturm-Liouville problem. A significant
 * feature of these polynomials is that they are orthogonal in the
 * interval \f$[-1,1]\f$ with respect to the function
 * \f$(1-x)^\alpha(1+x)^\beta (\alpha,\beta > -1)\f$
 *
 * Several functions related to the one-dimensional jacobi
 * polynomials: Evaluation, evaluation of derivatives, plus
 * computation of the roots via Newton's method.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see Karniadakis and Sherwin "Spectral/hp element methods for CFD"
 */
template<int N, typename T = double>
class Jacobi
{
public:

    /** @name Static values
     */
    //@{

    static const int order = N;
    static const int nOrder = N;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef Jacobi<N, T> self_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{


    /**
     * default values for a and b give the special case of Legendre
     * polynomials
     */
    Jacobi( value_type a = value_type( 0.0 ), value_type b = value_type( 0.0 ) )
        :
        M_a( a ),
        M_b( b )
    {}

    Jacobi( Jacobi const & p )
        :
        M_a( p.M_a ),
        M_b( p.M_b )
    {}

    ~Jacobi()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& p )
    {
        if ( this != &p )
        {
            M_a = p.M_a;
            M_b = p.M_b;
        }

        return *this;
    }

    /**
     * Evaluates the nth jacobi polynomial with weight parameters a,b
     * at a point x. Recurrence relations implemented from the
     * pseudocode given in Karniadakis and Sherwin, Appendix B
     *
     * a and b are defaulted to 0 and the Jacobi polynomial is then
     * the Legendre polynomial
     *
     * \param x point for polynomial evaluation
     * \return the value of the jacobi polynomial at \c x
     */
    value_type operator()( value_type const& x ) const;

    //@}

    /** @name Accessors
     */
    //@{

    uint16_type degree() const
    {
        return N;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Evaluates the nth jacobi polynomial with weight parameters a,b
     * at a point x. Recurrence relations implemented from the
     * pseudocode given in Karniadakis and Sherwin, Appendix B
     *
     * a and b are defaulted to 0 and the Jacobi polynomial is then
     * the Legendre polynomial
     *
     * \param x point for polynomial evaluation
     * \return the value of the jacobi polynomial at \c x
     */
    value_type value( value_type const& x ) const
    {
        return this->operator()( x );
    }

    value_type derivate( value_type const& x ) const;

    //@}



protected:

private:
    value_type M_a;
    value_type M_b;
};
template<int N, typename T>
typename Jacobi<N, T>::value_type
Jacobi<N, T>::operator()( value_type const& x ) const
{
    const value_type one = 1.0;
    const value_type two = 2.0;

    if ( N == 0 )
        return one;

    else if ( N == 1 )
        return 0.5 * ( M_a - M_b + ( M_a + M_b + two ) * x );

    else  // N >= 2
    {
        value_type apb = M_a + M_b;
        value_type pn2 = one;
        value_type pn1 = 0.5 * ( M_a - M_b + ( apb + two ) * x );
        value_type p = 0.0;

        for ( int k = 2; k < N+1; ++k )
        {
            value_type kv = value_type( k );
            value_type a1 = two * kv * ( kv + apb ) * ( two * kv + apb - two );
            value_type a2 = ( two * kv + apb - one ) * ( M_a * M_a - M_b * M_b );
            value_type a3 = ( ( two * kv + apb - two )
                              * ( two * kv + apb - one )
                              * ( two * kv + apb ) );
            value_type a4 = ( two * ( kv + M_a - one ) * ( kv + M_b - one )
                              * ( two * kv + apb ) );

            p = ( ( a2 + a3 * x ) * pn1 - a4 * pn2 )/a1;
            pn2 = pn1;
            pn1 = p;
        }

        return p;
    }
}
template<int N, typename T>
typename Jacobi<N, T>::value_type
Jacobi<N, T>::derivate( value_type const& x ) const
{
    if (  N == 0 )
        return 0.0;

    Jacobi<N-1, T> dp( M_a + 1.0, M_b + 1.0 );
    value_type Nv = value_type( N );
    return 0.5 * ( M_a  + M_b + Nv + 1.0 ) * dp( x );
}
namespace dyna
{
/**
 * \class Jacobi
 * \brief 1D Jacobi polynomial
 *
 * (excerpt from Karniadakis/Sherwin Appendix A)
 * Jacobi polynomials \f$ P^{\alpha,beta}_n(x)\f$ are a family of
 * polynomial to the singular Sturm-Liouville problem. A significant
 * feature of these polynomials is that they are orthogonal in the
 * interval \f$[-1,1]\f$ with respect to the function
 * \f$(1-x)^\alpha(1+x)^\beta (\alpha,\beta > -1)\f$
 *
 * Several functions related to the one-dimensional jacobi
 * polynomials: Evaluation, evaluation of derivatives, plus
 * computation of the roots via Newton's method.
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 * @see Karniadakis and Sherwin "Spectral/hp element methods for CFD"
 */
template<typename T = double>
class Jacobi
{
public:

    /** @name Static values
     */
    //@{

    //@}

    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef Jacobi<T> self_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{


    /**
     * default values for a and b give the special case of Legendre
     * polynomials
     */
    Jacobi( uint16_type N, value_type a = value_type( 0.0 ), value_type b = value_type( 0.0 ) )
        :
        M_degree( N ),
        M_a( a ),
        M_b( b )
    {}

    Jacobi( Jacobi const & p )
        :
        M_degree( p.M_degree ),
        M_a( p.M_a ),
        M_b( p.M_b )
    {}

    ~Jacobi()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type& operator=( self_type const& p )
    {
        if ( this != &p )
        {
            M_degree = p.M_degree;
            M_a = p.M_a;
            M_b = p.M_b;
        }

        return *this;
    }

    /**
     * Evaluates the nth jacobi polynomial with weight parameters a,b
     * at a point \p x. Recurrence relations implemented from the
     * pseudocode given in Karniadakis and Sherwin, Appendix B
     *
     * a and b are defaulted to 0 and the Jacobi polynomial is then
     * the Legendre polynomial
     *
     * \param x point for polynomial evaluation
     * \return the value of the jacobi polynomial at \c x
     */
    value_type operator()( value_type const& x ) const;

    //@}

    /** @name Accessors
     */
    //@{

    uint16_type degree() const
    {
        return M_degree;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    void setDegree( uint16_type N )
    {
        M_degree = N;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Evaluates the nth jacobi polynomial with weight parameters a,b
     * at a point x. Recurrence relations implemented from the
     * pseudocode given in Karniadakis and Sherwin, Appendix B
     *
     * a and b are defaulted to 0 and the Jacobi polynomial is then
     * the Legendre polynomial
     *
     * \param x point for polynomial evaluation
     * \return the value of the jacobi polynomial at \c x
     */
    value_type value( value_type const& x ) const
    {
        return this->operator()( x );
    }

    value_type derivate( value_type const& x ) const;

    //@}



protected:

private:
    uint16_type M_degree;
    value_type M_a;
    value_type M_b;
};
template<typename T>
typename Jacobi<T>::value_type
Jacobi<T>::operator()( value_type const& x ) const
{
    const uint16_type N = this->M_degree;
    const value_type one = 1.0;
    const value_type two = 2.0;

    if ( N == 0 )
        return one;

    else if ( N == 1 )
        return 0.5 * ( M_a - M_b + ( M_a + M_b + two ) * x );

    else  // N >= 2
    {
        value_type apb = M_a + M_b;
        value_type pn2 = one;
        value_type pn1 = 0.5 * ( M_a - M_b + ( apb + two ) * x );
        value_type p = 0.0;

        for ( uint16_type k = 2; k < N+1; ++k )
        {
            value_type kv = value_type( k );
            value_type a1 = two * kv * ( kv + apb ) * ( two * kv + apb - two );
            value_type a2 = ( two * kv + apb - one ) * ( M_a * M_a - M_b * M_b );
            value_type a3 = ( ( two * kv + apb - two )
                              * ( two * kv + apb - one )
                              * ( two * kv + apb ) );
            value_type a4 = ( two * ( kv + M_a - one ) * ( kv + M_b - one )
                              * ( two * kv + apb ) );

            p = ( ( a2 + a3 * x ) * pn1 - a4 * pn2 )/a1;
            pn2 = pn1;
            pn1 = p;
        }

        return p;
    }
}
template<typename T>
typename Jacobi<T>::value_type
Jacobi<T>::derivate( value_type const& x ) const
{
    const uint16_type N = this->M_degree;

    if (  N == 0 )
        return 0.0;

    Jacobi<T> dp( N-1, M_a + 1.0, M_b + 1.0 );
    value_type Nv = value_type( N );
    return 0.5 * ( M_a  + M_b + Nv + 1.0 ) * dp( x );
}
}

#if 0
template<uint16_type N, typename T>
typename Eigen::Matrix<T,N+1,Eigen::Dynamic>
JacobiBatchEvaluation( T a, T b, Matrix<T,Eigen::Dynamic,1> const& __pts )
{
    typedef T value_type;
    typedef Eigen::Matrix<T,Eigen::Dynamic,N+1> matrix_type;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> vector_type;
    matrix_type res( __pts.size(), N+1 );
    res.col( 0 ).setOnes();

    if ( N > 0 )
    {
        //ublas::col( res, 1 ) = 0.5 * ( ublas::scalar_vector<value_type>( res.size2(), a - b) + ( a + b + 2.0 ) * __pts );
        res.col( 1 ) = 0.5*( res.col( 1 ).Constant( res.rows(),a-b )+( a+b+2 )*__pts );
        value_type apb = a + b;

        for ( int k = 2; k < N+1; ++k )
        {
            value_type kv = value_type( k );
            value_type a1 = 2.0 * kv * ( kv + apb ) * ( 2.0 * kv + apb - 2.0 );
            value_type a2 = ( 2.0 * kv + apb - 1.0 ) * ( a * a - b * b );
            value_type a3 = ( 2.0 * kv + apb - 2.0 ) * ( 2.0 * kv + apb - 1.0 ) * ( 2.0 * kv + apb );
            value_type a4 = 2.0 * ( kv + a - 1.0 ) * ( kv + b - 1.0 ) * ( 2.0 * kv + apb );
            a2 = a2 / a1;
            a3 = a3 / a1;
            a4 = a4 / a1;

            res.col( k ) =
                ( a2* res.col( k-1 ) +
                  a3 * ( __pts.cwise()*res.col( k-1 ) ) -
                  a4 * res.col( k-2 ) );
        }
    }

    return res.transpose();
}
template<typename T,uint16_type N, uint16_type Np>
typename Eigen::Matrix<T,N+1,Np>
JacobiBatchEvaluation( T a, T b, Matrix<T,Np,1> const& __pts )
{
    typedef T value_type;
    typedef Eigen::Matrix<T,Np,N+1> matrix_type;
    typedef Eigen::Matrix<T,Np,1> vector_type;
    matrix_type res( Np, N+1 );
    res.col( 0 ).setOnes();

    if ( N > 0 )
    {
        //ublas::col( res, 1 ) = 0.5 * ( ublas::scalar_vector<value_type>( res.size2(), a - b) + ( a + b + 2.0 ) * __pts );
        res.col( 1 ) = 0.5*( res.col( 1 ).Constant( res.rows(),a-b )+( a+b+2 )*__pts );
        value_type apb = a + b;

        for ( int k = 2; k < N+1; ++k )
        {
            value_type kv = value_type( k );
            value_type a1 = 2.0 * kv * ( kv + apb ) * ( 2.0 * kv + apb - 2.0 );
            value_type a2 = ( 2.0 * kv + apb - 1.0 ) * ( a * a - b * b );
            value_type a3 = ( 2.0 * kv + apb - 2.0 ) * ( 2.0 * kv + apb - 1.0 ) * ( 2.0 * kv + apb );
            value_type a4 = 2.0 * ( kv + a - 1.0 ) * ( kv + b - 1.0 ) * ( 2.0 * kv + apb );
            a2 = a2 / a1;
            a3 = a3 / a1;
            a4 = a4 / a1;

            res.col( k ) =
                ( a2* res.col( k-1 ) +
                  a3 * ( __pts.cwise()*res.col( k-1 ) ) -
                  a4 * res.col( k-2 ) );
        }
    }

    return res.transpose();
}
#endif // 0 - Eigen

template<uint16_type N, typename T>
ublas::matrix<T>
JacobiBatchEvaluation( T a, T b, ublas::vector<T> const& __pts )
{
    typedef T value_type;
    ublas::matrix<T> res( N+1, __pts.size() );
    ublas::row( res, 0 ) = ublas::scalar_vector<value_type>( res.size2(), 1.0 );

    if ( N > 0 )
    {
        ublas::row( res, 1 ) = 0.5 * ( ublas::scalar_vector<value_type>( res.size2(), a - b ) + ( a + b + 2.0 ) * __pts );
        value_type apb = a + b;

        for ( int k = 2; k < N+1; ++k )
        {
            value_type kv = value_type( k );
            value_type a1 = 2.0 * kv * ( kv + apb ) * ( 2.0 * kv + apb - 2.0 );
            value_type a2 = ( 2.0 * kv + apb - 1.0 ) * ( a * a - b * b );
            value_type a3 = ( 2.0 * kv + apb - 2.0 ) * ( 2.0 * kv + apb - 1.0 ) * ( 2.0 * kv + apb );
            value_type a4 = 2.0 * ( kv + a - 1.0 ) * ( kv + b - 1.0 ) * ( 2.0 * kv + apb );
            a2 = a2 / a1;
            a3 = a3 / a1;
            a4 = a4 / a1;

            ublas::row( res, k ) = a2* ublas::row( res, k-1 ) +
                                   a3 * ublas::element_prod( __pts, ublas::row( res, k-1 ) ) -
                                   a4 * ublas::row( res, k-2 );
        }
    }

    return res;
}

template<uint16_type N, typename T>
void
JacobiBatchDerivation( ublas::matrix<T>& res, T a, T b, ublas::vector<T> const& __pts, mpl::bool_<true> )
{
    typedef T value_type;
    ublas::subrange( res, 1, N+1, 0, __pts.size() ) = JacobiBatchEvaluation<N-1, T>( a+1.0, b+1.0, __pts );

    for ( uint16_type i = 1; i < N+1; ++i )
        ublas::row( res, i ) *= 0.5*( a+b+value_type( i )+1.0 );
}

template<uint16_type N, typename T>
void
JacobiBatchDerivation( ublas::matrix<T>& /*res*/, T /*a*/, T /*b*/,
                       ublas::vector<T> const& /*__pts*/, mpl::bool_<false> )
{
}

template<uint16_type N, typename T>
ublas::matrix<T>
JacobiBatchDerivation( T a, T b, ublas::vector<T> const& __pts )
{
    typedef T value_type;
    ublas::matrix<T> res( N+1, __pts.size() );
    ublas::row( res, 0 ) = ublas::scalar_vector<value_type>( res.size2(), 0.0 );
    static const bool cond = N>0;
    JacobiBatchDerivation<N,T>( res, a, b, __pts, mpl::bool_<cond>() );
    return res;
}


#if 0
template<uint16_type N, typename T>
Eigen::Matrix<T,N+1,Eigen::Dynamic>
JacobiBatchDerivation( T a, T b, Eigen::Matrix<T,Eigen::Dynamic,1> const& __pts )
{
    typedef T value_type;
    Eigen::Matrix<T,N+1,Eigen::Dynamic> res( N+1, __pts.size() );
    res.row( 0 ).setZero();

    if ( N > 0 )
    {
        res.block( 1, 0, N, __pts.size() ) = JacobiBatchEvaluation<N-1, T>( a+1.0, b+1.0, __pts );

        //ublas::subrange( res, 1, N+1, 0, __pts.size() ) = JacobiBatchEvaluation<N-1, T>( a+1.0, b+1.0, __pts );
        for ( uint16_type i = 1; i < N+1; ++i )
            res.row( i ) *= 0.5*( a+b+value_type( i )+1.0 );
    }

    return res;
}
template<typename T,
         uint16_type N,
         uint16_type Ncols>
Eigen::Matrix<T,N+1,Ncols>
JacobiBatchDerivation( T a, T b, Eigen::Matrix<T,Ncols,1> const& __pts )
{
    typedef T value_type;
    Eigen::Matrix<T,N+1,Ncols> res;
    res.row( 0 ).setZero();

    if ( N > 0 )
    {
        res.block( 1, 0, N, __pts.size() ) = JacobiBatchEvaluation<T,N-1,Ncols>( a+1.0, b+1.0, __pts );

        //ublas::subrange( res, 1, N+1, 0, __pts.size() ) = JacobiBatchEvaluation<N-1, T>( a+1.0, b+1.0, __pts );
        for ( uint16_type i = 1; i < N+1; ++i )
            res.row( i ) *= 0.5*( a+b+value_type( i )+1.0 );
    }

    return res;
}
template<typename T>
Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
JacobiBatchDerivation( uint16_type N,
                       T a, T b, Eigen::Matrix<T,Eigen::Dynamic,1> const& __pts )
{
    typedef T value_type;
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> res( N+1, __pts.size() );
    res.row( 0 ).setZero();

    if ( N > 0 )
    {
        res.block( 1, 0, N, __pts.size() ) = JacobiBatchEvaluation<T>( N-1,
                                             a+1.0,
                                             b+1.0,
                                             __pts );

        for ( uint16_type i = 1; i < N+1; ++i )
            res.row( i ) *= 0.5*( a+b+value_type( i )+1.0 );
    }

    return res;
}
#endif // 0 - Eigen

namespace dyna
{
template<typename T>
ublas::matrix<T>
JacobiBatchEvaluation( uint16_type N, T a, T b, ublas::vector<T> const& __pts )
{
    typedef T value_type;
    ublas::matrix<T> res( N+1, __pts.size() );
    ublas::row( res, 0 ) = ublas::scalar_vector<value_type>( res.size2(), 1.0 );

    if ( N > 0 )
    {
        ublas::row( res, 1 ) = 0.5 * ( ublas::scalar_vector<value_type>( res.size2(), a - b ) + ( a + b + 2.0 ) * __pts );
        value_type apb = a + b;

        for ( uint16_type k = 2; k < N+1; ++k )
        {
            value_type kv = value_type( k );
            value_type a1 = 2.0 * kv * ( kv + apb ) * ( 2.0 * kv + apb - 2.0 );
            value_type a2 = ( 2.0 * kv + apb - 1.0 ) * ( a * a - b * b );
            value_type a3 = ( 2.0 * kv + apb - 2.0 ) * ( 2.0 * kv + apb - 1.0 ) * ( 2.0 * kv + apb );
            value_type a4 = 2.0 * ( kv + a - 1.0 ) * ( kv + b - 1.0 ) * ( 2.0 * kv + apb );
            a2 = a2 / a1;
            a3 = a3 / a1;
            a4 = a4 / a1;

            ublas::row( res, k ) = a2* ublas::row( res, k-1 ) +
                                   a3 * ublas::element_prod( __pts, ublas::row( res, k-1 ) ) -
                                   a4 * ublas::row( res, k-2 );
        }
    }

    return res;
}

template<typename T>
ublas::matrix<T>
JacobiBatchDerivation( uint16_type N, T a, T b, ublas::vector<T> const& __pts )
{
    typedef T value_type;
    ublas::matrix<T> res( N+1, __pts.size() );
    ublas::row( res, 0 ) = ublas::scalar_vector<value_type>( res.size2(), 0.0 );

    if ( N > 0 )
    {
        ublas::subrange( res, 1, N+1, 0, __pts.size() ) = JacobiBatchEvaluation<T>( N-1, a+1.0, b+1.0, __pts );

        for ( uint16_type i = 1; i < N+1; ++i )
            ublas::row( res, i ) *= 0.5*( a+b+value_type( i )+1.0 );
    }

    return res;
}
} // dyna

/// \cond DETAIL
namespace details
{

/**
 * Computes the m roots of \f$P_{m}^{a,b}\f$ on \f$[-1,1]\f$ by
 * Newton's method.  The initial guesses are the Chebyshev points as
 * we know an explicit formula for these polynomials.
 */
template<typename JacobiP, typename Vector>
void
roots( JacobiP const& p, Vector& xr )
{
    const uint16_type N = p.degree();
    typedef typename JacobiP::value_type value_type;

    if ( N != 0 )
    {
        value_type eps = type_traits<value_type>::epsilon();
        value_type r;
        int max_iter = 30;

        for ( int k = 0; k < N; ++k )
        {
            value_type pi = 4.0*math::atan( value_type( 1.0 ) );
            // use k-th checbychev point to  initiliaze newton
            r = -math::cos( ( 2.0*value_type( k ) + 1.0 ) * pi / ( 2.0 * value_type( N ) ) );

            // use average of r and xr[k-1] as starting point (see KS)
            if ( k > 0 )
                r = 0.5 * ( r + xr[k-1] );

            int j = 0;
            value_type jf = 2.0*eps;
            value_type delta  = 2.0*eps;

            do
            {
                // use deflation as proposed in KS
                value_type s = 0.0;

                for ( int i = 0; i < k; ++i )
                {
                    s +=  value_type( 1.0 ) / ( r - xr[i] );
                }

                jf = p( r );
                value_type jdf = p.derivate( r );
                delta = jf / ( jdf - jf * s );

                // newton step done
                r = r - delta;

                ++j;
            }
            while ( math::abs( jf ) > eps && j < max_iter );

            // store k-th root
            xr[k] = r;
        }
    }
}

template<typename T>
T fact( T  n )
{
    if ( n == 0.0 )
        return 1;

    return n*fact( n-1.0 );
}
template<int N, typename T, typename VectorW,  typename VectorN>
void
gaussjacobi( VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ) )
{
    typedef T value_type;

    Jacobi<N, T> p( a, b );

    roots( p, xr );

    const value_type two = 2.0;
    const value_type power = a+b+1.0;
    value_type a1 = math::pow( two,power );
    value_type a2 = fact( a+value_type( N ) );//gamma(a + m + 1);
    value_type a3 = fact( b+value_type( N ) );//gamma(b + m + 1);
    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact( value_type( N ) );
    value_type a6 = a1 * a2 * a3;// / a4 / a5;

    for ( int k = 0; k < N; ++k )
    {
        value_type fp = p.derivate( xr[k] );
        value_type dn = fp*fp*( 1.0  - xr[k]*xr[k] )*a4*a5;

        wr[k] =a6 / dn;
        //wr[k] = a6 / ( 1.0  - xr[k]*xr[k]) /  ( fp*fp );
    }
}


/**
 * \brief Gauss-Lobatto-Jacobi-Bouzitat Quadrature
 *
 * @see Karniadakis et Sherwin appendix B for theoretical details.
 **/
template<int N, typename T, typename VectorW,  typename VectorN>
void
gausslobattojacobi( VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ) )
{
    typedef T value_type;

    Jacobi<N-2, T> p( a+ 1.0, b+ 1.0 );

    Jacobi<N-1, T> q( a, b );

    VectorN prexr( N-2 );

    roots( p, prexr );

    xr[0]= -1.0;

    for ( int i = 1; i < N-1; ++i )
        xr[i] = prexr[i-1];

    xr[N-1]= 1.0;

    const value_type two = 2.0;

    value_type a1 = math::pow( two,int( a+b+1.0 ) );
    value_type a2 = fact( a+value_type( N ) -1.0 );//gamma(a + Q);
    value_type a3 = fact( b+value_type( N ) -1.0 );//gamma(b + Q);

    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact( value_type( N )-1.0 )*( value_type( N )-1.0 ); // (Q-1)!(Q-1)

    value_type a6 = a1 * a2 * a3;// Numerateur

    for ( int k = 0; k < N; ++k )
    {
        value_type fq = q.value( xr[k] );

        value_type dn = fq*fq*a4*a5;

        wr[k] =a6 /dn ;
    }

    wr[0] = wr[0]*( b+1.0 );
    wr[N-1] = wr[N-1]*( a+1.0 );
}

namespace dyna
{
/**
 * \brief Gauss-Lobatto-Jacobi-Bouzitat Quadrature
 *
 * @see Karniadakis et Sherwin appendix B for theoretical details.
 **/
template<typename T, typename VectorW,  typename VectorN>
void
gausslobattojacobi( size_type N, VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ), bool interior = false  )
{
    wr.resize( N - 2*interior );
    xr.resize( N - 2*interior );

    typedef T value_type;

    Feel::dyna::Jacobi<T> p( N-2, a+ 1.0, b+ 1.0 );

    Feel::dyna::Jacobi<T> q( N-1, a, b );

    VectorN prexr( N-2 );

    roots( p, prexr );

    if ( !interior )
    {
        xr[0]= -1.0;
        xr[N-1]= 1.0;
    }

    for ( size_type i = 1-interior; i < N-( 1+interior ); ++i )
        xr[i] = prexr[i-( 1-interior )];


    const value_type two = 2.0;

    value_type a1 = math::pow( two,int( a+b+1.0 ) );
    value_type a2 = fact( a+value_type( N ) -1.0 );//gamma(a + Q);
    value_type a3 = fact( b+value_type( N ) -1.0 );//gamma(b + Q);

    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact( value_type( N )-1.0 )*( value_type( N )-1.0 ); // (Q-1)!(Q-1)

    value_type a6 = a1 * a2 * a3;// Numerateur

    for ( int k = 1-interior; k < int( N-1-interior ); ++k )
    {
        value_type fq = q.value( xr[k] );

        value_type dn = fq*fq*a4*a5;

        wr[k] =a6 /dn ;
    }

    if ( !interior )
    {
        wr[0] = wr[0]*( b+1.0 );
        wr[N-1] = wr[N-1]*( a+1.0 );
    }
}

}

/** Implementation of Gauss-Radau-Jacobi-Bouzitat Quadrature
 *  You can see Karniadakis et Sherwin appendix B for theoretical details.
 *  I differentiate the choices of Gauss-Radau-Left and Gauss-Radau-Right
 *  quadrature methods.
 **/

template<int N, typename T, typename VectorW,  typename VectorN>
void
left_gaussradaujacobi( VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ) )
{
    typedef T value_type;

    Jacobi<N-1, T> p( a, b+1.0 );

    Jacobi<N-1, T> q( a, b );

    VectorN prexr( N-1 );

    roots( p, prexr );

    xr[0]= -1.0;

    for ( int i = 1; i < N; ++i )
        xr[i] = prexr[i-1];

    const value_type two = 2.0;

    value_type a1 = math::pow( two,int( a+b ) );
    value_type a2 = fact( a+value_type( N ) -1.0 );//gamma(a + Q);
    value_type a3 = fact( b+value_type( N ) -1.0 );//gamma(b + Q);

    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact( value_type( N )-1.0 )*( value_type( N )+b ); // (Q-1)!(Q+b)

    value_type a6 = a1 * a2 * a3;// Numerateur

    for ( int k = 0; k < N; ++k )
    {
        value_type fq = q.value( xr[k] );

        value_type dn = fq*fq*a4*a5;

        wr[k] =a6*( 1.0-xr[k] ) /dn ;
    }

    wr[0] = wr[0]*( b+1.0 );
}


template<int N, typename T, typename VectorW,  typename VectorN>
void
right_gaussradaujacobi( VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ) )
{
    typedef T value_type;

    Jacobi<N-1, T> p( a+1.0, b );

    Jacobi<N-1, T> q( a, b );

    VectorN prexr( N-1 );

    roots( p, prexr );

    for ( int i = 0; i < N-1; ++i )
        xr[i] = prexr[i];

    xr[N-1]=1.0;

    const value_type two = 2.0;

    value_type a1 = math::pow( two,int( a+b ) );
    value_type a2 = fact( a+value_type( N ) -1.0 );//gamma(a + Q);
    value_type a3 = fact( b+value_type( N ) -1.0 );//gamma(b + Q);

    value_type a4 = fact( a+b+value_type( N ) );//gamma(a + b + m + 1);
    value_type a5 = fact( value_type( N )-1.0 )*( value_type( N )+a ); // (Q-1)!(Q+a)

    value_type a6 = a1 * a2 * a3;// numerator

    for ( int k = 0; k < N; ++k )
    {
        value_type fq = q.value( xr[k] );

        value_type dn = fq*fq*a4*a5;

        wr[k] =a6*( 1.0+xr[k] ) /dn ;
    }

    wr[N-1] = wr[N-1]*( a+1.0 );
}


template<int N, typename T>
T integrate( boost::function<T( T const& )> const& f )
{
    typedef T value_type;
    ublas::vector<T> xr( Jacobi<N, T>::nOrder );
    ublas::vector<T> wr( Jacobi<N, T>::nOrder );

    // get weights and nodes for Legendre polynomials
    details::gaussjacobi<N, T, ublas::vector<T> >( wr, xr );

    value_type res = 0.0;

    for ( int k = 0; k < Jacobi<N, T>::nOrder; ++k )
        res += wr[k]*f( xr[k] );

    return res;
}
} // details
/// \endcond


} // Feel
#endif /* __Jacobi_H */
