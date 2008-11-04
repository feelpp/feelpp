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
/* File    :  gmm_precond_ildlt.h : modified version from I.T.L.           */
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
#ifndef GMM_PRECOND_ILDLT_H
#define GMM_PRECOND_ILDLT_H

// Incomplete Level 0 Cholesky Preconditioner.
// For use with symmetric real or hermitian complex sparse matrices.
//
// Notes: The idea under a concrete Preconditioner such 
//        as Incomplete Cholesky is to create a Preconditioner
//        object to use in iterative methods. 
//

// Y. Renard : Transformed in LDLT for stability reason.
//             U=LT is stored in csr format. D is stored on the diagonal of U.

#include <gmm_precond.h>

namespace gmm {

  template <typename Matrix>
  class ildlt_precond {

  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type magnitude_type;
    typedef csr_matrix_ref<value_type *, size_type *, size_type *, 0> tm_type;

    tm_type U;

  protected :
    std::vector<value_type> Tri_val;
    std::vector<size_type> Tri_ind, Tri_ptr;
 
    template<typename M> void do_ildlt(const M& A, row_major);
    void do_ildlt(const Matrix& A, col_major);

  public:

    size_type nrows(void) const { return mat_nrows(U); }
    size_type ncols(void) const { return mat_ncols(U); }
    value_type &D(size_type i) { return Tri_val[Tri_ptr[i]]; }
    const value_type &D(size_type i) const { return Tri_val[Tri_ptr[i]]; }
    ildlt_precond(void) {}
    void build_with(const Matrix& A) {
      Tri_ptr.resize(mat_nrows(A)+1);
      do_ildlt(A, typename principal_orientation_type<typename
		  linalg_traits<Matrix>::sub_orientation>::potype());
    }
    ildlt_precond(const Matrix& A)  { build_with(A); }
    size_type memsize() const { 
      return sizeof(*this) + 
	Tri_val.size() * sizeof(value_type) + 
	(Tri_ind.size()+Tri_ptr.size()) * sizeof(size_type); 
    }
  };

  template <typename Matrix> template<typename M>
  void ildlt_precond<Matrix>::do_ildlt(const M& A, row_major) {
    typedef typename linalg_traits<Matrix>::storage_type store_type;
    typedef value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    
    size_type Tri_loc = 0, n = mat_nrows(A), d, g, h, i, j, k;
    if (n == 0) return;
    T z, zz;
    Tri_ptr[0] = 0;
    R prec = default_tol(R());
    R max_pivot = gmm::abs(A(0,0)) * prec;
    
    for (int count = 0; count < 2; ++count) {
      if (count) { Tri_val.resize(Tri_loc); Tri_ind.resize(Tri_loc); }
      for (Tri_loc = 0, i = 0; i < n; ++i) {
	typedef typename linalg_traits<M>::const_sub_row_type row_type;
	row_type row = mat_const_row(A, i);
        typename linalg_traits<row_type>::const_iterator
	  it = vect_const_begin(row), ite = vect_const_end(row);

	if (count) { Tri_val[Tri_loc] = T(0); Tri_ind[Tri_loc] = i; }
	++Tri_loc; // diagonal element

	for (k = 0; it != ite; ++it, ++k) {
	  j = index_of_it(it, k, store_type());
	  if (i == j) {
	    if (count) Tri_val[Tri_loc-1] = *it; 
	  }
	  else if (j > i) {
	    if (count) { Tri_val[Tri_loc] = *it; Tri_ind[Tri_loc]=j; }
	    ++Tri_loc;
	  }
	}
	Tri_ptr[i+1] = Tri_loc;
      }
    }
    
    if (A(0,0) == T(0)) {
      Tri_val[Tri_ptr[0]] = T(1);
      DAL_WARNING(2, "pivot 0 is too small");
    }
    
    for (k = 0; k < n; k++) {
      d = Tri_ptr[k];
      z = T(gmm::real(Tri_val[d])); Tri_val[d] = z;
      if (gmm::abs(z) <= max_pivot) {
	Tri_val[d] = z = T(1);
	DAL_WARNING(2, "pivot " << k << " is too small");
      }
      max_pivot = std::max(max_pivot, std::min(gmm::abs(z) * prec, R(1)));
      
      for (i = d + 1; i < Tri_ptr[k+1]; ++i) Tri_val[i] /= z;
      for (i = d + 1; i < Tri_ptr[k+1]; ++i) {
	zz = gmm::conj(Tri_val[i] * z);
	h = Tri_ind[i];
	g = i;
	
	for (j = Tri_ptr[h] ; j < Tri_ptr[h+1]; ++j)
	  for ( ; g < Tri_ptr[k+1] && Tri_ind[g] <= Tri_ind[j]; ++g)
	    if (Tri_ind[g] == Tri_ind[j])
	      Tri_val[j] -= zz * Tri_val[g];
      }
    }
    U = tm_type(&(Tri_val[0]), &(Tri_ind[0]), &(Tri_ptr[0]),
			n, mat_ncols(A));
  }
  
  template <typename Matrix>
  void ildlt_precond<Matrix>::do_ildlt(const Matrix& A, col_major)
  { do_ildlt(gmm::conjugated(A), row_major()); }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const ildlt_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    gmm::lower_tri_solve(gmm::conjugated(P.U), v2, true);
    for (size_type i = 0; i < mat_nrows(P.U); ++i) v2[i] /= P.D(i);
    gmm::upper_tri_solve(P.U, v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const ildlt_precond<Matrix>& P,const V1 &v1,V2 &v2)
  { mult(P, v1, v2); }

  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const ildlt_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    gmm::lower_tri_solve(gmm::conjugated(P.U), v2, true);
    for (size_type i = 0; i < mat_nrows(P.U); ++i) v2[i] /= P.D(i);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const ildlt_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::upper_tri_solve(P.U, v2, true);  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const ildlt_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    copy(v1, v2);
    gmm::upper_tri_solve(P.U, v2, true);
    for (size_type i = 0; i < mat_nrows(P.U); ++i) v2[i] /= P.D(i);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const ildlt_precond<Matrix>& P, const V1 &v1,
			     V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::conjugated(P.U), v2, true); }



  // for compatibility with old versions

  template <typename Matrix>
  struct cholesky_precond : public ildlt_precond<Matrix> {
    cholesky_precond(const Matrix& A) : ildlt_precond<Matrix>(A) {}
    cholesky_precond(void) {}
  } IS_DEPRECATED;

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const cholesky_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    gmm::copy(v1, v2);
    gmm::lower_tri_solve(gmm::conjugated(P.U), v2, true);
    for (size_type i = 0; i < mat_nrows(P.U); ++i) v2[i] /= P.D(i);
    gmm::upper_tri_solve(P.U, v2, true);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const cholesky_precond<Matrix>& P,const V1 &v1,V2 &v2)
  { mult(P, v1, v2); }

  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const cholesky_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    copy(v1, v2);
    gmm::lower_tri_solve(gmm::conjugated(P.U), v2, true);
    for (size_type i = 0; i < mat_nrows(P.U); ++i) v2[i] /= P.D(i);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const cholesky_precond<Matrix>& P, const V1 &v1, V2 &v2)
  { copy(v1, v2); gmm::upper_tri_solve(P.U, v2, true);  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const cholesky_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    copy(v1, v2);
    gmm::upper_tri_solve(P.U, v2, true);
    for (size_type i = 0; i < mat_nrows(P.U); ++i) v2[i] /= P.D(i);
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const cholesky_precond<Matrix>& P, const V1 &v1,
			     V2 &v2)
  { copy(v1, v2); gmm::lower_tri_solve(gmm::conjugated(P.U), v2, true); }
  
}

#endif 

