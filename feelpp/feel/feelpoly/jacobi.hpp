/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2005-2006 EPFL
    SPDX-FileCopyrightText: 2006-2012 Joseph Fourier University
    SPDX-FileCopyrightText: 2012-2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
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

// clang-format off
#include <feel/feelcore/warnoff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <feel/feelcore/warnon.hpp>
// clang-format on

#include <feel/feelpoly/jacobirecurrence.hpp>
#include <feel/feelalg/eigenmap.hpp>

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
class Jacobi : public JacobiPolynomial<N, T>
{
    using super = JacobiPolynomial<N, T>;
public:
    /** Legacy compile-time degree constant. */
    static constexpr int order = N;
    /** Compile-time expansion degree. */
    static constexpr int nOrder = N;
    /** Scalar type retained for quadrature clients. */
    using value_type = T;
    /** Compatibility self type. */
    using self_type = Jacobi;
    using super::super;
    /** @brief Compatibility spelling of the mathematical derivative. */
    [[nodiscard]] T derivate( T const& x ) const { return this->derivative( x ); }
};

namespace dyna
{
/** @brief Compatibility facade for a runtime-degree Jacobi polynomial.
 * @tparam T Scalar type; the recurrence is shared with static polynomials.
 */
template<typename T = double>
class Jacobi : public JacobiPolynomial<Dynamic, T>
{
    using super = JacobiPolynomial<Dynamic, T>;
public:
    /** Scalar type retained for quadrature clients. */
    using value_type = T;
    /** Compatibility self type. */
    using self_type = Jacobi;
    /** @brief Construct a validated runtime-degree Jacobi polynomial. */
    Jacobi( int degree, T a = T( 0 ), T b = T( 0 ) )
        : super( RuntimeOrder{ degree }, a, b ) {}
    /** @brief Change the runtime degree after checking its supported range. */
    void setDegree( int degree ) { this->setOrderChecked( degree ); }
    /** @brief Compatibility spelling of the mathematical derivative. */
    [[nodiscard]] T derivate( T const& x ) const { return this->derivative( x ); }
};
} // namespace dyna

namespace detail
{
/** @brief Allocate a checked final legacy table before a shared kernel fills it. */
template<typename T>
[[nodiscard]] ublas::matrix<T> makeLegacyPolynomialTable( Eigen::Index rows, std::size_t samples )
{
    if ( samples > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
        throw std::length_error( "legacy polynomial sample count exceeds Eigen::Index" );
    (void)matrixElementCount( rows, static_cast<Eigen::Index>( samples ), sizeof( T ) );
    return ublas::matrix<T>( rows, samples );
}
} // namespace detail

/** @brief Legacy batch-value boundary forwarding to the Eigen recurrence. */
template<typename T>
ublas::matrix<T> JacobiBatchEvaluation( int degree, T a, T b, ublas::vector<T> const& points )
{
    JacobiPolynomial<Dynamic, T> polynomial{ RuntimeOrder{ degree }, a, b };
    auto result = ::Feel::detail::makeLegacyPolynomialTable<T>( polynomial.degree() + 1, points.size() );
    auto output = eigenMap( result );
    polynomial.valuesInto( output, eigenMap( points ) );
    return result;
}

/** @brief Legacy batch-derivative boundary forwarding to the shifted recurrence. */
template<typename T>
ublas::matrix<T> JacobiBatchDerivation( int degree, T a, T b, ublas::vector<T> const& points )
{
    JacobiPolynomial<Dynamic, T> polynomial{ RuntimeOrder{ degree }, a, b };
    auto result = ::Feel::detail::makeLegacyPolynomialTable<T>( polynomial.degree() + 1, points.size() );
    auto output = eigenMap( result );
    polynomial.derivativesInto( output, eigenMap( points ) );
    return result;
}

namespace dyna
{
/** @brief Runtime compatibility spelling of the shared value boundary. */
template<typename T>
ublas::matrix<T> JacobiBatchEvaluation( uint16_type degree, T a, T b, ublas::vector<T> const& points )
{ return Feel::JacobiBatchEvaluation( static_cast<int>( degree ), a, b, points ); }

/** @brief Runtime compatibility spelling of the shared derivative boundary. */
template<typename T>
ublas::matrix<T> JacobiBatchDerivation( uint16_type degree, T a, T b, ublas::vector<T> const& points )
{ return Feel::JacobiBatchDerivation( static_cast<int>( degree ), a, b, points ); }
} // namespace dyna

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
            // use k-th checbychev point to  initialize newton
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
template<typename T, typename VectorW,  typename VectorN>
void
gaussjacobi( int N, VectorW& wr, VectorN& xr, T a = T( 0.0 ), T b = T( 0.0 ) )
{
    typedef T value_type;

    Feel::dyna::Jacobi<T> p( N, a, b );

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
    Feel::details::gaussjacobi<N, T, ublas::vector<T> >( wr, xr );

    value_type res = 0.0;

    for ( int k = 0; k < Jacobi<N, T>::nOrder; ++k )
        res += wr[k]*f( xr[k] );

    return res;
}
} // details
/// \endcond


} // Feel
#endif /* __Jacobi_H */
