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
/* File    :  gmm_precond_ilu.h : modified version from I.T.L.             */
/*     									   */
/* Date : June 5, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
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
#ifndef GMM_PRECOND_ILU_H
#define GMM_PRECOND_ILU_H

//  Incomplete LU without fill-in Preconditioner.
//
// Notes: The idea under a concrete Preconditioner such 
//        as Incomplete LU is to create a Preconditioner
//        object to use in iterative methods. 
//

#include <gmm_precond.h>

namespace gmm {

  template <typename Matrix>
  class ilu_precond {

  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef csr_matrix_ref<value_type *, size_type *, size_type *, 0> tm_type;

    tm_type U, L;
    bool invert;
  protected :
    std::vector<value_type> L_val, U_val;
    std::vector<size_type> L_ind, U_ind, L_ptr, U_ptr;
 
    template<typename M> void do_ilu(const M& A, row_major);
    void do_ilu(const Matrix& A, col_major);

  public:
    
    size_type nrows(void) const { return mat_nrows(L); }
    size_type ncols(void) const { return mat_ncols(U); }
    
    void build_with(const Matrix& A) {
      invert = false;
       L_ptr.resize(mat_nrows(A)+1);
       U_ptr.resize(mat_nrows(A)+1);
       do_ilu(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
    ilu_precond(const Matrix& A) { build_with(A); }
    ilu_precond(void) {}
    size_type memsize() const { 
      return sizeof(*this) + 
	(L_val.size()+U_val.size()) * sizeof(value_type) + 
	(L_ind.size()+L_ptr.size()) * sizeof(size_type) +
	(U_ind.size()+U_ptr.size()) * sizeof(size_type); 
    }
  };

  template <typename Matrix> template <typename M>
  void ilu_precond<Matrix>::do_ilu(const M& A, row_major) {
    typedef typename linalg_traits<Matrix>::storage_type store_type;
    typedef value_type T;
    typedef typename number_traits<T>::magnitude_type R;

    size_type L_loc = 0, U_loc = 0, n = mat_nrows(A), i, j, k;
    if (n == 0) return;
    L_ptr[0] = 0; U_ptr[0] = 0;
    R prec = default_tol(R());
    R max_pivot = gmm::abs(A(0,0)) * prec;


    for (int count = 0; count < 2; ++count) {
      if (count) { 
	L_val.resize(L_loc); L_ind.resize(L_loc);
	U_val.resize(U_loc); U_ind.resize(U_loc);
      }
      L_loc = U_loc = 0;
      for (i = 0; i < n; ++i) {
	typedef typename linalg_traits<M>::const_sub_row_type row_type;
	row_type row = mat_const_row(A, i);
	typename linalg_traits<row_type>::const_iterator
	  it = vect_const_begin(row), ite = vect_const_end(row);
	
	if (count) { U_val[U_loc] = T(0); U_ind[U_loc] = i; }
	++U_loc; // diagonal element
	
	for (k = 0; it != ite; ++it, ++k) {
	  j = index_of_it(it, k, store_type());
	  if (j < i) {
	    if (count) { L_val[L_loc] = *it; L_ind[L_loc] = j; }
	    L_loc++;
	  }
	  else if (i == j) {
	    if (count) U_val[U_loc-1] = *it;
	  }
	  else {
	    if (count) { U_val[U_loc] = *it; U_ind[U_loc] = j; }
	    U_loc++;
	  }
	}
        L_ptr[i+1] = L_loc; U_ptr[i+1] = U_loc;
      }
    }
    
    if (A(0,0) == T(0)) {
      U_val[U_ptr[0]] = T(1);
      DAL_WARNING(2, "pivot 0 is too small");
    }

    size_type qn, pn, rn;
    for (i = 1; i < n; i++) {

      pn = U_ptr[i];
      if (gmm::abs(U_val[pn]) <= max_pivot) {
	U_val[pn] = T(1);
	DAL_WARNING(2, "pivot " << i << " is too small");
      }
      max_pivot = std::max(max_pivot,
			   std::min(gmm::abs(U_val[pn]) * prec, R(1)));

      for (j = L_ptr[i]; j < L_ptr[i+1]; j++) {
	pn = U_ptr[L_ind[j]];
	
	T multiplier = (L_val[j] /= U_val[pn]);
	
	qn = j + 1;
	rn = U_ptr[i];
	
	for (pn++; U_ind[pn] < i && pn < U_ptr[L_ind[j]+1]; pn++) {
	  while (L_ind[qn] < U_ind[pn] && qn < L_ptr[i+1])
	    qn++;
	  if (U_ind[pn] == L_ind[qn] && qn < L_ptr[i+1])
	    L_val[qn] -= multiplier * U_val[pn];
	}
	for (; pn < U_ptr[L_ind[j]+1]; pn++) {
	  while (U_ind[rn] < U_ind[pn] && rn < U_ptr[i+1])
	    rn++;
	  if (U_ind[pn] == U_ind[rn] && rn < U_ptr[i+1])
	    U_val[rn] -= multiplier * U_val[pn];
	}
      }
    }

    L = tm_type(&(L_val[0]), &(L_ind[0]), &(L_ptr[0]), n, mat_ncols(A));
    U = tm_type(&(U_val[0]), &(U_ind[0]), &(U_ptr[0]), n, mat_ncols(A));
  }
  
  template <typename Matrix>
  void ilu_precond<Matrix>::do_ilu(const Matrix& A, col_major) {
    do_ilu(gmm::transposed(A), row_major());
    invert = true;
  }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const ilu_precond<Matrix>& P, const V1 &v1, V2 &v2) {
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
  void transposed_mult(const ilu_precond<Matrix>& P,const V1 &v1,V2 &v2) {
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
  void left_mult(const ilu_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    else gmm::lower_tri_solve(P.L, v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const ilu_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    else gmm::upper_tri_solve(P.U, v2, false);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const ilu_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::upper_tri_solve(P.U, v2, false);
    else gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const ilu_precond<Matrix>& P, const V1 &v1,
			     V2 &v2) {
    copy(v1, v2);
    if (P.invert) gmm::lower_tri_solve(P.L, v2, true);
    else gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
  }


}

#endif 

