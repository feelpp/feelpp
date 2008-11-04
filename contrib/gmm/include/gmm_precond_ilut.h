/* -*- c++ -*- (enables emacs c++ mode)                                    */
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
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_ilut.h : modified version from I.T.L.            */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@insa-toulouse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003-2004  Yves Renard.                                   */
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
#ifndef GMM_PRECOND_ILUT_H
#define GMM_PRECOND_ILUT_H

//: ILUT:  Incomplete LU with threshold and K fill-in Preconditioner.
//  The algorithm of ILUT(A, 0, 1.0e-6) is slower than ILU(A). If No fill-in 
//  is arrowed, you can use ILU instead of ILUT.
//
// Notes: The idea under a concrete Preconditioner such 
//        as ilut is to create a Preconditioner
//        object to use in iterative methods. 
//

/*
  Performane comparing for SSOR, ILU and ILUT based on sherman 5 matrix 
  in Harwell-Boeing collection on Sun Ultra 30 UPA/PCI (UltraSPARC-II 296MHz)
  Preconditioner & Factorization time  &  Number of Iteration \\ \hline
  SSOR        &   0.010577  & 41 \\
  ILU         &   0.019336  & 32 \\
  ILUT with 0 fill-in and threshold of 1.0e-6 & 0.343612 &  23 \\
  ILUT with 5 fill-in and threshold of 1.0e-6 & 0.343612 &  18 \\ \hline
*/

#include <gmm_precond.h>

namespace gmm {

  template<typename T> struct elt_rsvector_value_less_ {
    inline bool operator()(const elt_rsvector_<T>& a, 
			   const elt_rsvector_<T>& b) const
    { return (gmm::abs(a.e) > gmm::abs(b.e)); }
  };

  template <typename Matrix>
  class ilut_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef rsvector<value_type> svector;
    typedef row_matrix<svector> LU_Matrix;

    bool invert;
    LU_Matrix L, U;

  protected:
    int K;
    double eps;    

    template<typename M> void do_ilut(const M&, row_major);
    void do_ilut(const Matrix&, col_major);

  public:
    void build_with(const Matrix& A) {
      invert = false;
      gmm::resize(L, mat_nrows(A), mat_ncols(A));
      gmm::resize(U, mat_nrows(A), mat_ncols(A));
      do_ilut(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
    ilut_precond(const Matrix& A, int k_, double eps_) 
      : L(mat_nrows(A), mat_ncols(A)), U(mat_nrows(A), mat_ncols(A)),
	K(k_), eps(eps_) { build_with(A); }
    ilut_precond(int k_, double eps_) :  K(k_), eps(eps_) {}
    ilut_precond(void) { K = 10; eps = 1E-7; }
    size_type memsize() const { 
      return sizeof(*this) + (nnz(U)+nnz(L))*sizeof(value_type);
    }
  };

  template<typename Matrix> template<typename M> 
  void ilut_precond<Matrix>::do_ilut(const M& A, row_major) {
    typedef value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    
    size_type n = mat_nrows(A);
    if (n == 0) return;
    std::vector<T> indiag(n);
    svector w(mat_ncols(A)), wL(mat_ncols(A)), wU(mat_ncols(A));
    T tmp;
    gmm::clear(U); gmm::clear(L);
    R prec = default_tol(R()); 
    R max_pivot = gmm::abs(A(0,0)) * prec;

    for (size_type i = 0; i < n; ++i) {
      gmm::copy(mat_const_row(A, i), w);
      double norm_row = gmm::vect_norm2(w);

      size_type nL = 0, nU = 1;
      if (is_sparse(A)) {
	typename linalg_traits<svector>::iterator it = vect_begin(w),
	  ite = vect_end(w);
	for (; it != ite; ++it) if (i > it.index()) nL++;
	nU = w.nb_stored() - nL;
      }

      for (size_type krow = 0, k; krow < w.nb_stored(); ++krow) {
	typename svector::iterator wk = w.begin() + krow;
	if ((k = wk->c) >= i) break;
	tmp = (wk->e) * indiag[k];
	if (gmm::abs(tmp) < eps * norm_row) { w.sup(k); --krow; } 
	else { wk->e += tmp; gmm::add(scaled(mat_row(U, k), -tmp), w); }
      }
      tmp = w[i];

      if (gmm::abs(tmp) <= max_pivot) {
	DAL_WARNING(2, "pivot " << i << " too small. try with ilutp ?");
	w[i] = tmp = T(1);
      }

      max_pivot = std::max(max_pivot, std::min(gmm::abs(tmp) * prec, R(1)));
      indiag[i] = T(1) / tmp;
      gmm::clean(w, eps * norm_row);
      std::sort(w.begin(), w.end(), elt_rsvector_value_less_<T>());
      typename svector::const_iterator wit = w.begin(), wite = w.end();
      size_type nnl = 0, nnu = 0;
      
      wL.base_resize(nL+K); wU.base_resize(nU+K+1);
      typename svector::iterator witL = wL.begin(), witU = wU.begin();
      for (; wit != wite; ++wit) 
	if (wit->c < i) { if (nnl < nL+K) { *witL++ = *wit; ++nnl; } }
	else { if (nnu < nU+K  || wit->c == i) { *witU++ = *wit; ++nnu; } }
      wL.base_resize(nnl); wU.base_resize(nnu);
      std::sort(wL.begin(), wL.end());
      std::sort(wU.begin(), wU.end());
      gmm::copy(wL, L.row(i));
      gmm::copy(wU, U.row(i));
      
//       wit = w.begin(); nnl = 0; nnu = 0;
//       for (; wit != wite; ++wit) // copy to be optimized ...
//   	if (wit->c < i) { if (nnl < nL+K) { L(i, wit->c) = wit->e; ++nnl;} }
//   	else if (nnu < nU+K || wit->c == i) { U(i, wit->c) = wit->e; ++nnu; }
    }

  }

  template<typename Matrix> 
  void ilut_precond<Matrix>::do_ilut(const Matrix& A, col_major) {
    do_ilut(gmm::transposed(A), row_major());
    invert = true;
  }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const ilut_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    if (P.invert) {
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
    else {
      gmm::lower_tri_solve(P.L, v2, true);
      gmm::upper_tri_solve(P.U, v2, false);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const ilut_precond<Matrix>& P,const V1 &v1,V2 &v2) {
    gmm::copy(v1, v2);
    if (P.invert) {
      gmm::lower_tri_solve(P.L, v2, true);
      gmm::upper_tri_solve(P.U, v2, false);
    }
    else {
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const ilut_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    else gmm::lower_tri_solve(P.L, v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const ilut_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    else gmm::upper_tri_solve(P.U, v2, false);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const ilut_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::upper_tri_solve(P.U, v2, false);
    else gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const ilut_precond<Matrix>& P, const V1 &v1,
			     V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::lower_tri_solve(P.L, v2, true);
    else gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
  }

}

#endif 

