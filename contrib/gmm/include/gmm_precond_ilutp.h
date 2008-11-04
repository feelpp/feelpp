/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_precond_ilutp.h : ILUTP preconditionner for sparse       */
/*                                  matrices                               */
/*     									   */
/* Date : October 14, 2004.                                                */
/* Author : Yves Renard, Yves.Renard@insa-toulouse.fr                      */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2004  Yves Renard.                                        */
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
#ifndef GMM_PRECOND_ILUTP_H
#define GMM_PRECOND_ILUTP_H

// ILUTP:  Incomplete LU with threshold and K fill-in Preconditioner and
//         column pivoting (See Yousef Saad, Iterative Methods for
//         sparse linear systems, PWS Publishing Company, section 10.4.4

// TODO : store the permutation by cycles to avoid the temporary vector

#include <gmm_precond_ilut.h>

namespace gmm {

  template <typename Matrix>
  class ilutp_precond  {
  public :
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef rsvector<value_type> svector;
    typedef row_matrix<svector> LU_Matrix;
    typedef col_matrix<svector> CLU_Matrix;

    bool invert;
    LU_Matrix L, U;
    gmm::unsorted_sub_index indperm;
    gmm::unsorted_sub_index indperminv;    
    mutable std::vector<value_type> temporary;

  protected:
    int K;
    double eps;

    template<typename M> void do_ilutp(const M&, row_major);
    void do_ilutp(const Matrix&, col_major);

  public:
    void build_with(const Matrix& A) {
      invert = false;
      gmm::resize(L, mat_nrows(A), mat_ncols(A));
      gmm::resize(U, mat_nrows(A), mat_ncols(A));
      do_ilutp(A, typename principal_orientation_type<typename
	      linalg_traits<Matrix>::sub_orientation>::potype());
    }
    ilutp_precond(const Matrix& A, int k_, double eps_) 
      : L(mat_nrows(A), mat_ncols(A)), U(mat_nrows(A), mat_ncols(A)),
	K(k_), eps(eps_) { build_with(A); }
    ilutp_precond(int k_, double eps_) :  K(k_), eps(eps_) {}
    ilutp_precond(void) { K = 10; eps = 1E-7; }
    size_type memsize() const { 
      return sizeof(*this) + (nnz(U)+nnz(L))*sizeof(value_type);
    }
  };


  template<typename Matrix> template<typename M> 
  void ilutp_precond<Matrix>::do_ilutp(const M& A, row_major) {
    typedef value_type T;
    typedef typename number_traits<T>::magnitude_type R;
    
    size_type n = mat_nrows(A);
    CLU_Matrix CU(n,n);
    if (n == 0) return;
    std::vector<T> indiag(n);
    temporary.resize(n);
    std::vector<size_type> ipvt(n), ipvtinv(n);
    for (size_type i = 0; i < n; ++i) ipvt[i] = ipvtinv[i] = i;
    indperm = unsorted_sub_index(ipvt);
    indperminv = unsorted_sub_index(ipvtinv);
    svector w(mat_ncols(A));
    
    T tmp = T(0);
    gmm::clear(L); gmm::clear(U);
    R prec = default_tol(R()); 
    R max_pivot = gmm::abs(A(0,0)) * prec;

    for (size_type i = 0; i < n; ++i) {
      
      copy(sub_vector(mat_const_row(A, i), indperm), w);
      double norm_row = gmm::vect_norm2(mat_const_row(A, i)); 

      size_type nL = 0, nU = 1;
      if (is_sparse(A))
	{ nL = nnz(mat_const_row(A, i)) / 2; nU = nL + 1; }

      for (size_type krow = 0, k; krow < w.nb_stored(); ++krow) {
	typename svector::iterator wk = w.begin() + krow;
	if ((k = wk->c) >= i) break;
	tmp = (wk->e) * indiag[k];
	if (gmm::abs(tmp) < eps * norm_row) { w.sup(k); --krow; } 
	else { wk->e += tmp; gmm::add(scaled(mat_row(U, k), -tmp), w); }
      }
      gmm::clean(w, eps * norm_row);
      std::sort(w.begin(), w.end(), elt_rsvector_value_less_<T>());
      typename svector::const_iterator wit = w.begin(), wite = w.end();
      size_type ip = size_type(-1);
      for (; wit != wite; ++wit)
	if (wit->c >= i) { ip = wit->c; tmp = wit->e; break; }
      if (ip == size_type(-1) || gmm::abs(tmp) <= max_pivot)
	{ DAL_WARNING(2, "pivot " << i << " too small"); ip=i; w[i]=tmp=T(1); }
      max_pivot = std::max(max_pivot, std::min(gmm::abs(tmp) * prec, R(1)));
      indiag[i] = T(1) / tmp;
      wit = w.begin();
      size_type nnl = 0, nnu = 0;
      for (; wit != wite; ++wit) {
	if (wit->c < i) { if (nnl < nL+K) { L(i, wit->c) = wit->e; ++nnl; } }
	else if (nnu < nU+K) { CU(i, wit->c) = U(i, wit->c) = wit->e; ++nnu; }
      }
      if (ip != i) {
	typename svector::const_iterator iti = CU.col(i).begin();
	typename svector::const_iterator itie = CU.col(i).end();
	typename svector::const_iterator itp = CU.col(ip).begin();
	typename svector::const_iterator itpe = CU.col(ip).end();
	
	while (iti != itie && itp != itpe) {
	  if (iti->c < itp->c) { U.row(iti->c).swap_indices(i, ip); ++iti; }
	  else if (iti->c > itp->c) { U.row(itp->c).swap_indices(i,ip);++itp; }
	  else { U.row(iti->c).swap_indices(i, ip); ++iti; ++itp; }
	}
	for( ; iti != itie; ++iti) U.row(iti->c).swap_indices(i, ip);
	for( ; itp != itpe; ++itp) U.row(itp->c).swap_indices(i, ip);

	CU.swap_col(i, ip);
	
	indperm.swap(i, ip);
	indperminv.swap(ipvt[i], ipvt[ip]);
	std::swap(ipvtinv[ipvt[i]], ipvtinv[ipvt[ip]]);
	std::swap(ipvt[i], ipvt[ip]);
      }
    }
  }

  template<typename Matrix> 
  void ilutp_precond<Matrix>::do_ilutp(const Matrix& A, col_major) {
    do_ilutp(gmm::transposed(A), row_major());
    invert = true;
  }

  template <typename Matrix, typename V1, typename V2> inline
  void mult(const ilutp_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.invert) {
      gmm::copy(gmm::sub_vector(v1, P.indperm), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
    else {
      gmm::copy(v1, P.temporary);
      gmm::lower_tri_solve(P.L, P.temporary, true);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperminv), v2);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_mult(const ilutp_precond<Matrix>& P,const V1 &v1,V2 &v2) {
    if (P.invert) {
      gmm::copy(v1, P.temporary);
      gmm::lower_tri_solve(P.L, P.temporary, true);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperminv), v2);
    }
    else {
      gmm::copy(gmm::sub_vector(v1, P.indperm), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void left_mult(const ilutp_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.invert) {
      gmm::copy(gmm::sub_vector(v1, P.indperm), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    }
    else {
      copy(v1, v2);
      gmm::lower_tri_solve(P.L, v2, true);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void right_mult(const ilutp_precond<Matrix>& P, const V1 &v1, V2 &v2) {
    if (P.invert) {
      copy(v1, v2);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
    else {
      copy(v1, P.temporary);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperminv), v2);
    }
  }

  template <typename Matrix, typename V1, typename V2> inline
  void transposed_left_mult(const ilutp_precond<Matrix>& P, const V1 &v1,
			    V2 &v2) {
    if (P.invert) {
      copy(v1, P.temporary);
      gmm::upper_tri_solve(P.U, P.temporary, false);
      gmm::copy(gmm::sub_vector(P.temporary, P.indperminv), v2);
    }
    else {
      copy(v1, v2);
      gmm::upper_tri_solve(gmm::transposed(P.L), v2, true);
    }
  }
  
  template <typename Matrix, typename V1, typename V2> inline
  void transposed_right_mult(const ilutp_precond<Matrix>& P, const V1 &v1,
			     V2 &v2) {
    if (P.invert) {
      copy(v1, v2);
      gmm::lower_tri_solve(P.L, v2, true);
    }
    else {
      gmm::copy(gmm::sub_vector(v1, P.indperm), v2);
      gmm::lower_tri_solve(gmm::transposed(P.U), v2, false);
    }
  }

}

#endif 

