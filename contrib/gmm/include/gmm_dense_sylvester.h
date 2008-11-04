/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_dense_sylvester.h : solve sylvester equation for dense   */
/*                                    systems.                             */
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

#ifndef GMM_DENSE_SYLVESTER_H
#define GMM_DENSE_SYLVESTER_H

#include <gmm_kernel.h>

namespace gmm {



  /* ********************************************************************* */
  /*   Kronecker system matrix.                                            */
  /* ********************************************************************* */
  template <typename MAT1, typename MAT2, typename MAT3>
  void kron(const MAT1 &m1, const MAT2 &m2, const MAT3 &m3_,
	    bool init = true) {
    MAT3 &m3 = const_cast<MAT3 &>(m3_);
    size_type m = mat_nrows(m1), n = mat_ncols(m1);
    size_type l = mat_nrows(m2), k = mat_ncols(m2);

    if (mat_nrows(m3) != m*l || mat_ncols(m3) != n*k)
      DAL_THROW(dimension_error, "dimensions mismatch");

    for (size_type i = 0; i < m; ++i)
      for (size_type j = 0; j < m; ++j)
	if (init)
	  gmm::copy(gmm::scaled(m2, m1(i,j)),
		    gmm::sub_matrix(m3, sub_interval(l*i, l),
				    sub_interval(k*j, k)));
	else
	  gmm::add(gmm::scaled(m2, m1(i,j)),
		    gmm::sub_matrix(m3, sub_interval(l*i, l),
				    sub_interval(k*j, k)));
  }
	

  /* ********************************************************************* */
  /*   Copy a matrix into a vector.                                        */
  /* ********************************************************************* */

  template <typename MAT, typename VECT>
  colmatrix_to_vector(const MAT &A, VECT &v, col_major) {
    size_type m = mat_nrows(A), n = mat_ncols(A);
    if (m*n != vect_size(v)) DAL_THROW(dimension_error, "dimensions mismatch");
    for (size_type i = 0; i < n; ++i)
      gmm::copy(mat_col(A, i), sub_vector(v, sub_interval(i*m, m)));
  }

  template <typename MAT, typename VECT>
  colmatrix_to_vector(const MAT &A, VECT &v, row_and_col)
  { colmatrix_to_vector(A, v, col_major()); }

  template <typename MAT, typename VECT>
  colmatrix_to_vector(const MAT &A, VECT &v, col_and_row)
  { colmatrix_to_vector(A, v, col_major()); }

  template <typename MAT, typename VECT>
  colmatrix_to_vector(const MAT &A, VECT &v, row_major) {
    size_type m = mat_nrows(mat), n = mat_ncols(A);
    if (m*n != vect_size(v)) DAL_THROW(dimension_error, "dimensions mismatch");
    for (size_type i = 0; i < m; ++i)
      gmm::copy(mat_row(A, i), sub_vector(v, sub_slice(i, n, m)));
  }

  template <typename MAT, typename VECT> inline
  colmatrix_to_vector(const MAT &A, const VECT &v_) {
    VECT &v = const_cast<VECT &>(v_);
    colmatrix_to_vector(A, v, typename linalg_traits<MAT>::sub_orientation());
  }


  /* ********************************************************************* */
  /*   Copy a vector into a matrix.                                        */
  /* ********************************************************************* */

  template <typename MAT, typename VECT>
  vector_to_colmatrix(const VECT &v, MAT &A, col_major) {
    size_type m = mat_nrows(A), n = mat_ncols(A);
    if (m*n != vect_size(v)) DAL_THROW(dimension_error, "dimensions mismatch");
    for (size_type i = 0; i < n; ++i)
      gmm::copy(sub_vector(v, sub_interval(i*m, m)), mat_col(A, i));
  }

  template <typename MAT, typename VECT>
  vector_to_colmatrix(const VECT &v, MAT &A, row_and_col)
  { vector_to_colmatrix(v, A, col_major()); }

  template <typename MAT, typename VECT>
  vector_to_colmatrix(const VECT &v, MAT &A, col_and_row)
  { vector_to_colmatrix(v, A, col_major()); }

  template <typename MAT, typename VECT>
  vector_to_colmatrix(const VECT &v, MAT &A, row_major) {
    size_type m = mat_nrows(mat), n = mat_ncols(A);
    if (m*n != vect_size(v)) DAL_THROW(dimension_error, "dimensions mismatch");
    for (size_type i = 0; i < m; ++i)
      gmm::copy(sub_vector(v, sub_slice(i, n, m)), mat_row(A, i));
  }

  template <typename MAT, typename VECT> inline
  vector_to_colmatrix(const VECT &v, const MAT &A_) {
    MAT &A = const_cast<MAT &>(A_);
    vector_to_colmatrix(v, A, typename linalg_traits<MAT>::sub_orientation());
  }

  /* ********************************************************************* */
  /*   Solve sylvester equation.                                           */
  /* ********************************************************************* */

  // very prohibitive solver, to be replaced ... 
  template <typename MAT1, typename MAT2, typename MAT3, typename MAT4 >
  void sylvester(const MAT1 &m1, const MAT2 &m2, const MAT3 &m3,
		 const MAT4 &m4_) {
    typedef typename linalg_traits<Mat>::value_type T;
    
    MAT3 &m4 = const_cast<MAT4 &>(m4_);
    size_type m = mat_nrows(m1), n = mat_ncols(m1);
    size_type l = mat_nrows(m2), k = mat_ncols(m2);
    
    if (m != n || l != k || m != mat_nrows(m3)
	|| l != mat_ncols(m3) || m != mat_nrows(m4) || l != mat_ncols(m4))
      DAL_THROW(dimension_error, "dimensions mismatch");

    gmm::dense_matrix<T> akronb(m*l, m*l);
    gmm::dense_matrix<T> idm(m, m), idl(l,l);
    gmm::copy(identity_matrix(), idm);
    gmm::copy(identity_matrix(), idl);
    std::vector<T> x(m*l), c(m*l);
    
    kron(idl, m1, akronb);
    kron(gmm::transposed(m2), idm, akronb, false);

    colmatrix_to_vector(m3, c);
    lu_solve(akronb, c, x);
    vector_to_colmatrix(x, m4);

  }
}

#endif

