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
/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solver_gmres.h : from I.T.L.                             */
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

#ifndef GMM_KRYLOV_GMRES_H
#define GMM_KRYLOV_GMRES_H

#include <gmm_kernel.h>
#include <gmm_iter.h>
#include <gmm_modified_gram_schmidt.h>

namespace gmm {

  // Generalized Minimum Residual
  //
  //   This solve the unsymmetric linear system Ax = b using restarted GMRES.
  //
  //   See: Y. Saad and M. Schulter. GMRES: A generalized minimum residual
  //   algorithm for solving nonsysmmetric linear systems, SIAM
  //   J. Sci. Statist. Comp.  7(1986), pp, 856-869
  //

  template <typename Mat, typename Vec, typename VecB, typename Precond,
	    typename Basis >
  void gmres(const Mat &A, Vec &x, const VecB &b, const Precond &M,
	     int restart, iteration &outer, Basis& KS) {

    typedef typename linalg_traits<Vec>::value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    std::vector<T> w(vect_size(x)), r(vect_size(x)), u(vect_size(x));
    std::vector<T> c_rot(restart+1), s_rot(restart+1), s(restart+1);
    gmm::dense_matrix<T> H(restart+1, restart);

    mult(M,b,r);
    outer.set_rhsnorm(gmm::vect_norm2(r));
    if (outer.get_rhsnorm() == 0.0) { clear(x); return; }
    
    mult(A, scaled(x, -1), b, w);
    mult(M, w, r);
    R beta = gmm::vect_norm2(r), beta_old = beta;
    int blocked = 0;

    iteration inner = outer;
    inner.reduce_noisy();
    inner.set_maxiter(restart);
    inner.set_name("GMRes inner");

    while (! outer.finished(beta)) {
      
      gmm::copy(gmm::scaled(r, R(1)/beta), KS[0]);
      gmm::clear(s);
      s[0] = beta;
      
      size_type i = 0; inner.init();
      
      do {
	mult(A, KS[i], u);
	mult(M, u, KS[i+1]);
	orthogonalize(KS, mat_col(H, i), i);
	R a = gmm::vect_norm2(KS[i+1]);
	H(i+1, i) = T(a);
	gmm::scale(KS[i+1], T(1) / a);
	for (size_type k = 0; k < i; ++k)
	  Apply_Givens_rotation_left(H(k,i), H(k+1,i), c_rot[k], s_rot[k]);
	
	Givens_rotation(H(i,i), H(i+1,i), c_rot[i], s_rot[i]);
	Apply_Givens_rotation_left(H(i,i), H(i+1,i), c_rot[i], s_rot[i]);
	Apply_Givens_rotation_left(s[i], s[i+1], c_rot[i], s_rot[i]);
	
	++inner, ++outer, ++i;
      } while (! inner.finished(gmm::abs(s[i])));

      upper_tri_solve(H, s, i, false);
      combine(KS, s, x, i);
      mult(A, gmm::scaled(x, -1), b, w);
      mult(M, w, r);
      beta_old = std::min(beta, beta_old); beta = gmm::vect_norm2(r);
      if (int(inner.get_iteration()) < restart -1 || beta_old <= beta)
	++blocked; else blocked = 0;
      if (blocked > 10) {
	if (outer.get_noisy()) cout << "Gmres is blocked, exiting\n";
	break;
      }
    }
  }


  template <typename Mat, typename Vec, typename VecB, typename Precond >
  void gmres(const Mat &A, Vec &x, const VecB &b,
	     const Precond &M, int restart, iteration& outer) {
    typedef typename linalg_traits<Vec>::value_type T;
    modified_gram_schmidt<T> orth(restart, vect_size(x));
    gmres(A, x, b, M, restart, outer, orth); 
  }

}

#endif
