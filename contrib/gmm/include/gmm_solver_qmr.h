// -*- c++ -*-
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu>
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// You should have received a copy of the License Agreement for the
// Iterative Template Library along with the software;  see the
// file LICENSE.
//
// Permission to modify the code and to distribute modified code is
// granted, provided the text of this NOTICE is retained, a notice that
// the code was modified is included with the above COPYRIGHT NOTICE and
// with the COPYRIGHT NOTICE in the LICENSE file, and that the LICENSE
// file is distributed with the modified code.
//
// LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
// By way of example, but not limitation, Licensor MAKES NO
// REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
// PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
// OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
// OR OTHER RIGHTS.
//=======================================================================
//
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solver_qmr.h : from I.T.L.                               */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GMM++                                            */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */
#ifndef GMM_QMR_H
#define GMM_QMR_H

#include <gmm_kernel.h>
#include <gmm_iter.h>

namespace gmm {

  //: Quasi-Minimal Residual
  //
  //  This routine solves the unsymmetric linear system Ax = b using the
  //  Quasi-Minimal Residual method.
  //
  //See: R. W. Freund and N. M. Nachtigal, A quasi-minimal residual method for
  //non-Hermitian linear systems, Numerical Math., 60(1991), pp. 315-339
  //
  // Preconditioner -  Incomplete LU, Incomplete LU with threshold,
  //                   SSOR or identity_preconditioner.

  template <typename Matrix, typename Vector, typename VectorB,
	    typename Precond1>
  void qmr(const Matrix &A, Vector &x, const VectorB &b, const Precond1 &M1,
	   iteration& iter) {

    typedef typename linalg_traits<Vector>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    T delta(0), ep(0), beta(0), theta_1(0), gamma_1(0);
    T theta(0), gamma(1), eta(-1);
    R rho_1(0), rho, xi;

    typedef typename temporary_vector<Vector>::vector_type TmpVec;
    size_type nn = vect_size(x);
    TmpVec r(nn), v_tld(nn), y(nn), w_tld(nn), z(nn), v(nn), w(nn);
    TmpVec y_tld(nn), z_tld(nn), p(nn), q(nn), p_tld(nn), d(nn), s(nn);

    iter.set_rhsnorm(double(gmm::vect_norm2(b)));
    if (iter.get_rhsnorm() == 0.0) { clear(x); return; }

    gmm::mult(A, gmm::scaled(x, T(-1)), b, r);
    gmm::copy(r, v_tld);

    gmm::left_mult(M1, v_tld, y);
    rho = gmm::vect_norm2(y);

    gmm::copy(r, w_tld);
    gmm::transposed_right_mult(M1, w_tld, z);
    xi = gmm::vect_norm2(z);

    while (! iter.finished_vect(r)) {

        if (rho == R(0) || xi == R(0))
            {
                if (iter.get_maxiter() == size_type(-1))
                    { DAL_THROW(failure_error, "QMR failed to converge"); }
                else { DAL_WARNING(1, "QMR failed to converge"); return; }
            }
      gmm::copy(gmm::scaled(v_tld, T(R(1)/rho)), v);
      gmm::scale(y, T(R(1)/rho));

      gmm::copy(gmm::scaled(w_tld, T(R(1)/xi)), w);
      gmm::scale(z, T(R(1)/xi));

      delta = gmm::vect_sp(z, y);
      if (delta == T(0))
          {
              if (iter.get_maxiter() == size_type(-1))
                  { DAL_THROW(failure_error, "QMR failed to converge"); }
              else { DAL_WARNING(1, "QMR failed to converge"); return; }
          }
      gmm::right_mult(M1, y, y_tld);
      gmm::transposed_left_mult(M1, z, z_tld);

      if (iter.first()) {
	gmm::copy(y_tld, p);
	gmm::copy(z_tld, q);
      } else {
	gmm::add(y_tld, gmm::scaled(p, -(T(xi  * delta) / ep)), p);
	gmm::add(z_tld, gmm::scaled(q, -(T(rho * delta) / ep)), q);
      }

      gmm::mult(A, p, p_tld);

      ep = gmm::vect_sp(q, p_tld);
      if (ep == T(0))
          {
              if (iter.get_maxiter() == size_type(-1))
                  { DAL_THROW(failure_error, "QMR failed to converge"); }
              else { DAL_WARNING(1, "QMR failed to converge"); return; }
          }

      beta = ep / delta;
      if (beta == T(0))
          {
              if (iter.get_maxiter() == size_type(-1))
                  { DAL_THROW(failure_error, "QMR failed to converge"); }
              else { DAL_WARNING(1, "QMR failed to converge"); return; }
          }

      gmm::add(p_tld, gmm::scaled(v, -beta), v_tld);
      gmm::left_mult(M1, v_tld, y);

      rho_1 = rho;
      rho = gmm::vect_norm2(y);

      gmm::mult(gmm::transposed(A), q, w_tld);
      gmm::add(w_tld, gmm::scaled(w, -beta), w_tld);
      gmm::transposed_right_mult(M1, w_tld, z);

      xi = gmm::vect_norm2(z);

      gamma_1 = gamma;
      theta_1 = theta;

      theta = rho / (gamma_1 * beta);
      gamma = T(1) / gmm::sqrt(T(1) + gmm::sqr(theta));

      if (gamma == T(0))
          {
              if (iter.get_maxiter() == size_type(-1))
                  { DAL_THROW(failure_error, "QMR failed to converge"); }
              else { DAL_WARNING(1, "QMR failed to converge"); return; }
          }

      eta = -eta * T(rho_1) * gmm::sqr(gamma) / (beta * gmm::sqr(gamma_1));

      if (iter.first()) {
	gmm::copy(gmm::scaled(p, eta), d);
	gmm::copy(gmm::scaled(p_tld, eta), s);
      } else {
	T tmp = gmm::sqr(theta_1 * gamma);
	gmm::add(gmm::scaled(p, eta), gmm::scaled(d, tmp), d);
	gmm::add(gmm::scaled(p_tld, eta), gmm::scaled(s, tmp), s);
      }
      gmm::add(d, x);
      gmm::add(gmm::scaled(s, T(-1)), r);

      ++iter;
    }
  }


}

#endif

