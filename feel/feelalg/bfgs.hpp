/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-05-02

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
   \file bfgs.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-05-02
 */
#include <feel/feelcore/typetraits.hpp>

#ifndef _BFGS_HPP_
#define _BFGS_HPP_

namespace Feel
{
namespace ublas = boost::numeric::ublas;
/**
 * BFGS algorithm (Broyden, Fletcher, Goldfarb, Shanno)
 * Quasi Newton method for optimization problems.
 * with Wolfe Line search.
 *
 * Ripped from getfem++ by Y. Renard
 */
enum BFGSType { BFGS = 0,  DFP };

/**
   delta[k] = x[k+1] - x[k]<br>
   gamma[k] = grad f(x[k+1]) - grad f(x[k])<br>
   H[0] = I<br>
   BFGS : zeta[k] = delta[k] - H[k] gamma[k]<br>
   DFP  : zeta[k] = H[k] gamma[k]<br>
   tau[k] = gamma[k]^T zeta[k]<br>
   rho[k] = 1 / gamma[k]^T delta[k]<br>
   BFGS : H[k+1] = H[k] + rho[k](zeta[k] delta[k]^T + delta[k] zeta[k]^T)<br>
   - rho[k]^2 tau[k] delta[k] delta[k]^T<br>
   DFP  : H[k+1] = H[k] + rho[k] delta[k] delta[k]^T<br>
   - (1/tau[k])zeta[k] zeta[k]^T
*/
// Object representing the inverse of the Hessian
template <typename VECTOR>
struct BFGSInvHessian
{
    typedef VECTOR vector_type;
    typedef typename vector_type::value_type T;
    typedef typename vector_type::value_type value_type;
    //typedef typename number_traits<T>::magnitude_type R;
    typedef value_type magnitude_type;


    BFGSInvHessian( BFGSType v = BFGS )
    {
        version = v;
    }

    template<typename VEC1, typename VEC2>
    void hmult( const VEC1 &X, VEC2 &Y )
    {
        Y.assign(  X );

        for ( size_type k = 0 ; k < delta.size(); ++k )
        {
            T xdelta = ublas::inner_prod( X, delta[k] );
            T xzeta = ublas::inner_prod( X, zeta[k] );

            switch ( version )
            {
            case BFGS :
                Y.plus_assign( rho[k]*xdelta*zeta[k] );
                Y.plus_assign( rho[k]*( xzeta-rho[k]*tau[k]*xdelta )*delta[k] );
                break;

            case DFP :
                Y.plus_assign( rho[k]*xdelta*delta[k] );
                Y.minus_assign( xzeta/tau[k]*zeta[k] );
                break;
            }
        }
    }

    void restart( void )
    {
        delta.clear();
        gamma.clear();
        zeta.clear();
        tau.clear();
        rho.clear();
    }

    template<typename VECT1, typename VECT2>
    void update( const VECT1 &deltak, const VECT2 &gammak )
    {
        size_type N = deltak.size();
        size_type k = delta.size();

        vector_type Y( N );

        hmult( gammak, Y );

        delta.resize( k+1 );
        gamma.resize( k+1 );
        zeta.resize( k+1 );
        tau.resize( k+1 );
        rho.resize( k+1 );

        delta[k].resize( N );
        gamma[k].resize( N );
        zeta[k].resize( N );

        delta[k].assign( deltak );
        gamma[k].assign( gammak );

        rho[k] = magnitude_type( 1 ) / ublas::inner_prod( deltak, gammak );

        if ( version == BFGS )
            zeta[k].plus_assign( delta[k] - Y );

        else
            zeta[k].assign( Y );

        tau[k] = ublas::inner_prod( gammak,  zeta[k] );
    }

    //
    // Data
    //
    std::vector<vector_type> delta, gamma, zeta;
    std::vector<T> tau, rho;
    int version;
};


template <typename FUNCTION, typename DERIVATIVE, typename VECTOR,  typename IterationBFGS>
void bfgs( FUNCTION f,
           DERIVATIVE grad,
           VECTOR &x,
           int restart,
           IterationBFGS& iter,
           BFGSType version = BFGS,
           float lambda_init = 0.001,
           float /*print_norm*/ = 1.0 )
{

    //typedef typename linalg_traits<VECTOR>::value_type T;
    //typedef typename number_traits<T>::magnitude_type R;

    typedef VECTOR vector_type;
    typedef typename vector_type::value_type T;
    typedef typename vector_type::value_type value_type;
    typedef typename type_traits<T>::real_type real_type;
    typedef value_type magnitude_type;


    BFGSInvHessian<vector_type> invhessian( version );

    VECTOR r( x.size() ), d( x.size() ), y( x.size() ), r2( x.size() );
    grad( x, r );
    real_type lambda = lambda_init, valx = f( x ), valy;
    int nb_restart( 0 );

    //if (iter.get_noisy() >= 1) cout << "value " << valx / print_norm << " ";
    while ( ! iter.isFinished( r ) )
    {
        invhessian.hmult( r, d );
        d *= -1;

        // Wolfe Line search
        real_type derivative = ublas::inner_prod( r, d );
        real_type lambda_min( 0 );
        real_type lambda_max( 0 );
        real_type m1 = 0.27;
        real_type m2 = 0.57;
        bool unbounded = true, blocked = false, grad_computed = false;

        for ( ;; )
        {

            //add(x, scaled(d, lambda), y);
            y = lambda*d+x;
            valy = f( y );

#if 0

            if ( iter.get_noisy() >= 2 )
            {
                cout << "Wolfe line search, lambda = " << lambda
                     << " value = " << valy /print_norm << endl;
            }

#endif

            if ( valy <= valx + m1 * lambda * derivative )
            {
                grad( y, r2 );
                grad_computed = true;

                T derivative2 = ublas::inner_prod( r2, d );

                if ( derivative2 >= m2*derivative )
                    break;

                lambda_min = lambda;
            }

            else
            {
                lambda_max = lambda;
                unbounded = false;
            }

            if ( unbounded )
                lambda *= real_type( 10 );

            else
                lambda = ( lambda_max + lambda_min ) / real_type( 2 );

            if ( valy <= real_type( 2 )*valx &&
                    ( lambda < real_type( lambda_init*1E-8 ) ||
                      ( !unbounded && lambda_max-lambda_min < real_type( lambda_init*1E-8 ) ) ) )
            {
                blocked = true;
                lambda = lambda_init;
                break;
            }
        }

        // Rank two update
        ++iter;

        if ( !grad_computed )
            grad( y, r2 );

        //gmm::add(scaled(r2, -1), r);
        r.minus_assign( r2 );

        if ( iter.numberOfIterations() % restart == 0 || blocked )
        {
            //if (iter.get_noisy() >= 1) cout << "Restart\n";
            invhessian.restart();

            if ( ++nb_restart > 10 )
            {
                //if (iter.get_noisy() >= 1) cout << "BFGS is blocked, exiting\n";
                return;
            }
        }

        else
        {
            //invhessian.update(gmm::scaled(d,lambda), gmm::scaled(r,-1));
            invhessian.update( lambda*d, -r );
            nb_restart = 0;
        }

        r.assign( r2 );
        x.assign( y );
        valx = valy;
        //if (iter.get_noisy() >= 1) cout << "BFGS value " << valx/print_norm << "\t";
    }

}


template <typename FUNCTION, typename DERIVATIVE, typename VECTOR,  typename IterationBFGS>
inline void
dfp( FUNCTION f,
     DERIVATIVE grad,
     VECTOR &x,
     int restart,
     IterationBFGS& iter,
     BFGSType version = DFP )
{
    bfgs( f, grad, x, restart, iter, version );
}


}
#endif
