/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 12 Sep 2014

   Copyright (C) 2014 Feel++ Consortium

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/*
  this work has been adapted from Shogun
 */
#ifndef FEELPP_JACOBI_ELLIPTIC_FUNCTIONS_HPP
#define FEELPP_JACOBI_ELLIPTIC_FUNCTIONS_HPP

#include <limits>
#include <cmath>

#include <feel/feelcore/feel.hpp>


namespace Feel
{

/** @brief Class that contains methods for computing Jacobi elliptic functions
 * related to complex analysis. These functions are inverse of the elliptic
 * integral of first kind, i.e.
 * \f[
 * u(k,m)=\int_{0}^{k}\frac{dt}{\sqrt{(1-t^{2})(1-m^{2}t^{2})}}
 * =\int_{0}^{\varphi}\frac{d\theta}{\sqrt{(1-m^{2}sin^{2}\theta)}}
 * \f]
 * where \f$k=sin\varphi\f$, \f$t=sin\theta\f$ and parameter \f$m, 0\le m
 * \le 1\f$ is called modulus. Three main Jacobi elliptic functions are defined
 * as \f$sn(u,m)=k=sin\theta\f$, \f$cn(u,m)=cos\theta=\sqrt{1-sn(u,m)^{2}}\f$
 * and \f$dn(u,m)=\sqrt{1-m^{2}sn(u,m)^{2}}\f$.
 * For \f$k=1\f$, i.e. \f$\varphi=\frac{\pi}{2}\f$, \f$u(1,m)=K(m)\f$ is known
 * as the complete elliptic integral of first kind. Similarly, \f$u(1,m'))=
 * K'(m')\f$, \f$m'=\sqrt{1-m^{2}}\f$ is called the complementary complete
 * elliptic integral of first kind. Jacobi functions are double periodic with
 * quardratic periods \f$K\f$ and \f$K'\f$.
 *
 * This namespace provides two sets of methods for computing \f$K,K'\f$, and
 * \f$sn,cn,dn\f$. Useful for computing rational approximation of matrix
 * functions given by Cauchy's integral formula, etc.
 */
namespace math
{
    namespace detail
    {
        template<typename T>
        inline T computeQuarterPeriod(T b)
        {
            const T eps=std::numeric_limits<T>::epsilon();
            const T pi=M_PI;

            T a=1.0;
            T mm=1.0;

            int64_t p=2;
            do
            {
                T a_new=(a+b)*0.5;
                T b_new=sqrt(a*b);
                T c=(a-b)*0.5;
                mm=T(p)*c*c;
                p<<=1;
                a=a_new;
                b=b_new;
            } while (mm>eps);
            return pi*0.5/a;
        }

        template<typename T>
        inline T polySix(T x)
        {
            return (132*pow(x,6)+42*pow(x,5)+14*pow(x,4)+5*pow(x,3)+2*pow(x,2)+x);
        }
    }

	/** Computes the quarter periods (K and K') of Jacobian elliptic functions
	 * (see class description).
	 * @param L
	 * @param K the quarter period (to be computed) on the T axis
	 * @param Kp the quarter period (to be computed) on the Imaginary axis
	 * computed
	 */
    template<typename T>
    void ellipkkp(T L, T &K, T &Kp)
    {
        CHECK(L>=0.0) <<"ellipKKp(): Parameter L should be non-negative\n";
        const T eps=std::numeric_limits<T>::epsilon();
        const T pi=M_PI;

        if (L>10.0)
        {
            K=pi*0.5;
            Kp=pi*L+log(4.0);
        }
        else
        {
            T m=exp(-2.0*pi*L);
            T mp=1.0-m;
            if (m<eps)
            {
                K=detail::computeQuarterPeriod(sqrt(mp));
                Kp=T(std::numeric_limits<float64_t>::max());
            }
            else if (mp<eps)
            {
                K=T(std::numeric_limits<float64_t>::max());
                Kp=detail::computeQuarterPeriod(sqrt(m));
            }
            else
            {
                K=detail::computeQuarterPeriod(sqrt(mp));
                Kp=detail::computeQuarterPeriod(sqrt(m));
            }
        }
    }

	/** Computes three main Jacobi elliptic functions, \f$sn(u,m)\f$,
	 * \f$cn(u,m)\f$ and \f$dn(u,m)\f$ (see class description).
	 * @param u the elliptic integral of the first kind \f$u(k,m)\f$
	 * @param m the modulus parameter, \f$0\le m \le 1\f$
	 * @param sn Jacobi elliptic function sn(u,m)
	 * @param cn Jacobi elliptic function cn(u,m)
	 * @param dn Jacobi elliptic function dn(u,m)
	 */
    template<typename T>
    void ellipjc(std::complex<T> u, T m, std::complex<T> &sn, std::complex<T> &cn,
                 std::complex<T> &dn)
    {
        CHECK(m>=0.0 && m<=1.0) <<		"ellipjc(): \
		Parameter m should be >=0 and <=1\n";

        const T eps=sqrt(std::numeric_limits<T>::epsilon());

        if (m>=(1.0-eps))
        {
            std::complex<T> t=std::tanh(u);
            std::complex<T> b=std::cosh(u);
            std::complex<T> ai=0.25*(1.0-m);
            std::complex<T> twon=b*std::sinh(u);
            sn=t+ai*(twon-u)/(b*b);
            std::complex<T> phi=T(1.0)/b;
            ai*=t*phi;
            cn=phi-ai*(twon-u);
            dn=phi+ai*(twon+u);
        }
        else
        {
            const T prec=4.0*eps;
            const int32_type MAX_ITER=128;
            int32_type i=0;
            T kappa[MAX_ITER];

            while (i<MAX_ITER && m>prec)
            {
                T k;
                if (m>0.001)
                {
                    T mp=sqrt(1.0-m);
                    k=(1.0-mp)/(1.0+mp);
                }
                else
                    k=detail::polySix(m/4.0);
                u/=(1.0+k);
                m=k*k;
                kappa[i++]=k;
            }
            std::complex<T> sin_u=sin(u);
            std::complex<T> cos_u=cos(u);
            std::complex<T> t=T(0.25*m)*(u-sin_u*cos_u);
            sn=sin_u-t*cos_u;
            cn=cos_u+t*sin_u;
            dn=T(1.0)+T(0.5*m)*(cos_u*cos_u);

            i--;
            while (i>=0)
            {
                T k=kappa[i--];
                std::complex<T> ksn2=k*(sn*sn);
                std::complex<T> d=T(1.0)+ksn2;
                sn*=(1.0+k)/d;
                cn*=dn/d;
                dn=(T(1.0)-ksn2)/d;
            }
        }
    }

    template<typename T>
    void ellipjc(std::vector<std::complex<T>> u, T l, std::vector<std::complex<T>> &sn,
                 std::vector<std::complex<T>> &cn, std::vector<std::complex<T>> &dn, bool flag=true)
     {

         const T eps=sqrt(std::numeric_limits<T>::epsilon());
         std::vector<size_type> high;
         double m;

         if ( flag )
         {
             double k, kp;
             math::ellipkkp(l, k, kp);

             for (size_type i = 0; i < u.size(); ++i)
             {
                 if (std::imag(u[i]) > kp/2)
                     high.push_back(i);
             }

             for (size_type const& s : high )
             {
                 u[s] = std::complex<T>(0,kp) - u[s];
             }

             m = std::exp(-2*pi*l);
         }
         else
         {
             m = l;
         }

         if ( m < 4*eps )
         {
             for (size_type i = 0; i < u.size(); ++i)
             {
                 sn.push_back(std::sin(u[i]) + (m/4)*(std::sin(u[i])*std::cos(u[i])-u[i])*std::cos(u[i]));
                 cn.push_back(std::cos(u[i]) + (m/4)*(-std::sin(u[i])*std::cos(u[i])+u[i])*std::sin(u[i]));
                 dn.push_back(std::complex<T>(1,0) + (m/4)*(std::pow(std::cos(u[i]),2.)-std::pow(std::sin(u[i]),2.)-std::complex<T>(1,0)));
             }
         }
         else
         {
             double kappa = 0;

             if ( m > 1e-03 )
             {
                 kappa = (1-sqrt(1-m))/(1+sqrt(1-m));
             }
             else
             {
                 std::vector<double> polycoeffs = {132, 42, 14, 5, 2, 1, 0};
                 for ( size_type expt=0; expt<polycoeffs.size(); ++expt )
                 {
                     kappa += polycoeffs[expt]*std::pow((m/4),polycoeffs.size()-expt-1);
                 }
             }

             double mu = std::pow(kappa,2.);
             std::vector<std::complex<T>> v;
             for (size_type i = 0; i < u.size(); ++i)
             {
                 v.push_back((1/(1+kappa))*u[i]);
             }

             std::vector<std::complex<T>> sn1, cn1, dn1;
             ellipjc(v, mu, sn1, cn1, dn1, false);

             std::vector<std::complex<T>> denom;
             for (size_type i = 0; i < u.size(); ++i)
             {
                 denom.push_back(std::complex<T>(1,0)+kappa*std::pow(sn1[i],2.));
             }

             for (size_type i = 0; i < u.size(); ++i)
             {
                 sn.push_back((1+kappa)*sn1[i]/denom[i]);
                 cn.push_back(cn1[i]*dn1[i]/denom[i]);
                 dn.push_back((std::complex<T>(1,0)-kappa*std::pow(sn1[i],2))/denom[i]);
             }

             if (!high.empty())
             {
                 std::vector<std::complex<T>> snh(sn.size());
                 std::vector<std::complex<T>> cnh(sn.size());
                 std::vector<std::complex<T>> dnh(sn.size());

                 snh = sn;
                 cnh = cn;
                 dnh = dn;

                 for (size_type const& s : high )
                 {
                     sn[s] = -std::complex<T>(1,0)/(sqrt(m)*sn[s]);
                     cn[s] = std::complex<T>(0,1)*dnh[s]/(sqrt(m)*sn[s]);
                     dn[s] = std::complex<T>(0,1)*cnh[s]/(sn[s]);
                 }
             }
         }

     }
}

}

#endif /* FEELPP_JACOBI_ELLIPTIC_FUNCTIONS_HPP */
