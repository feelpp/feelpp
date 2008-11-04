/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library : Generic Matrix Methods  (gmm)                                 */
/* File    : gmm_solver_constrained_cg.h : constrained conjugate, gradient */
/*           Modified version of I.T.L. conjugated gradient.               */
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
//=======================================================================
// Copyright (C) 1997-2001
// Authors: Andrew Lumsdaine <lums@osl.iu.edu> 
//          Lie-Quan Lee     <llee@osl.iu.edu>
//
// This file is part of the Iterative Template Library
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


//  preconditionning does not work

#ifndef GMM_SOLVER_CCG_H__
#define GMM_SOLVER_CCG_H__

#include <gmm_kernel.h>
#include <gmm_iter.h>

namespace gmm {

  template <typename CMatrix, typename CINVMatrix, typename Matps,
	    typename VectorX>
  void pseudo_inverse(const CMatrix &C, CINVMatrix &CINV,
		      const Matps& /* PS */, VectorX&) {
    // compute the pseudo inverse of the non-square matrix C such
    // CINV = inv(C * trans(C)) * C.
    // based on a conjugate gradient method.
    
    // optimisable : copie de la ligne, precalcul de C * trans(C).
    
    typedef VectorX TmpVec;
    typedef size_t size_type;
    typedef typename linalg_traits<VectorX>::value_type value_type;
    
    size_type nr = mat_nrows(C), nc = mat_ncols(C);
    
    TmpVec d(nr), e(nr), l(nc), p(nr), q(nr), r(nr);
    value_type rho, rho_1, alpha;
    clear(d);
    clear(CINV);
    
    for (size_type i = 0; i < nr; ++i) {
      d[i] = 1.0; rho = 1.0;
      clear(e);
      copy(d, r);
      copy(d, p);
      
      while (rho >= 1E-38) { /* conjugate gradient to compute e             */
	                     /* which is the i nd row of inv(C * trans(C))  */
	mult(gmm::transposed(C), p, l);
	mult(C, l, q);	  
	alpha = rho / vect_sp(p, q);
	add(scaled(p, alpha), e);  
	add(scaled(q, -alpha), r); 
	rho_1 = rho;
	rho = vect_sp(r, r);
	add(r, scaled(p, rho / rho_1), p);
      }
      
      mult(transposed(C), e, l); /* l is the i nd row of CINV     */
      // cout << "l = " << l << endl;
      clean(l, 1E-15);
      copy(l, mat_row(CINV, i));
      
      d[i] = 0.0;
    }
  }
  
  // Compute the minimum of 1/2((Ax).x) - bx under the contraint Cx <= f

  template < typename Matrix,  typename CMatrix, typename Matps,
	     typename VectorX, typename VectorB, typename VectorF,
	     typename Preconditioner >
  void constrained_cg(const Matrix& A, const CMatrix& C, VectorX& x,
		      const VectorB& b, const VectorF& f,const Matps& PS,
		      const Preconditioner& M, iteration &iter) {
    typedef typename temporary_dense_vector<VectorX>::vector_type TmpVec;
    typedef typename temporary_vector<CMatrix>::vector_type TmpCVec;
    typedef row_matrix<TmpCVec> TmpCmat;
    
    typedef size_t size_type;
    typedef typename linalg_traits<VectorX>::value_type value_type;
    value_type rho = 1.0, rho_1, lambda, gamma;
    TmpVec p(vect_size(x)), q(vect_size(x)), q2(vect_size(x)),
      r(vect_size(x)), old_z(vect_size(x)), z(vect_size(x)),
      memox(vect_size(x));
    std::vector<bool> satured(mat_nrows(C));
    clear(p);
    iter.set_rhsnorm(sqrt(vect_sp(PS, b, b)));
    if (iter.get_rhsnorm() == 0.0) iter.set_rhsnorm(1.0);
   
    TmpCmat CINV(mat_nrows(C), mat_ncols(C));
    pseudo_inverse(C, CINV, PS, x);
    
    while(true) {
      // computation of residu
      copy(z, old_z);
      copy(x, memox);
      mult(A, scaled(x, -1.0), b, r);
      mult(M, r, z); // preconditionner not coherent
      bool transition = false;
      for (size_type i = 0; i < mat_nrows(C); ++i) {
	value_type al = vect_sp(mat_row(C, i), x) - f[i];
	if (al >= -1.0E-15) {
	  if (!satured[i]) { satured[i] = true; transition = true; }
	  value_type bb = vect_sp(mat_row(CINV, i), z);
	  if (bb > 0.0) add(scaled(mat_row(C, i), -bb), z);
	}
	else
	  satured[i] = false;
      }
    
      // descent direction
      rho_1 = rho; rho = vect_sp(PS, r, z); // ...
      
      if (iter.finished(rho)) break;
      
      if (iter.get_noisy() > 0 && transition) std::cout << "transition\n";
      if (transition || iter.first()) gamma = 0.0;
      else gamma = std::max(0.0, (rho - vect_sp(PS, old_z, z) ) / rho_1);
      // std::cout << "gamma = " << gamma << endl;
      // itl::add(r, itl::scaled(p, gamma), p);
      add(z, scaled(p, gamma), p); // ...
      
      ++iter;
      // one dimensionnal optimization
      mult(A, p, q);
      lambda = rho / vect_sp(PS, q, p);
      for (size_type i = 0; i < mat_nrows(C); ++i)
	if (!satured[i]) {
	  value_type bb = vect_sp(mat_row(C, i), p) - f[i];
	  if (bb > 0.0)
	    lambda = std::min(lambda, (f[i]-vect_sp(mat_row(C, i), x)) / bb);
	}
      add(x, scaled(p, lambda), x);
      add(memox, scaled(x, -1.0), memox);
      
    }
  }
  
}

#endif //  GMM_SOLVER_CCG_H__
