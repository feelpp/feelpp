/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_solvers_Schwarz_additive.h : generic solver.             */
/*     									   */
/* Date : October 13, 2002.                                                */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr                      */
/*          Michel Fournie, fournie@mip.ups-tlse.fr                        */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2002-2004  Yves Renard, Michel Fournié.                   */
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


#ifndef GMM_SOLVERS_SCHWARZ_ADDITIVE_H__
#define GMM_SOLVERS_SCHWARZ_ADDITIVE_H__

#include <gmm_kernel.h>
#include <gmm_superlu_interface.h>
#include <gmm_solver_cg.h>
#include <gmm_solver_gmres.h>
#include <gmm_solver_bicgstab.h>
#include <gmm_solver_qmr.h>

namespace gmm {
      
  /* ******************************************************************** */
  /*		Additive Schwarz interfaced local solvers                 */
  /* ******************************************************************** */

  struct using_cg {};
  struct using_gmres {};
  struct using_bicgstab {};
  struct using_qmr {};

  template <typename P, typename local_solver, typename Matrix>
  struct actual_precond {
    typedef P APrecond;
    static const APrecond &transform(const P &PP) { return PP; }
  };

  template <typename Matrix1, typename Precond, typename Vector> 
  void AS_local_solve(using_cg, const Matrix1 &A, Vector &x, const Vector &b,
		 const Precond &P, iteration &iter)
  { cg(A, x, b, P, iter); }

  template <typename Matrix1, typename Precond, typename Vector> 
  void AS_local_solve(using_gmres, const Matrix1 &A, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { gmres(A, x, b, P, 100, iter); }
  
  template <typename Matrix1, typename Precond, typename Vector> 
  void AS_local_solve(using_bicgstab, const Matrix1 &A, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { bicgstab(A, x, b, P, iter); }

  template <typename Matrix1, typename Precond, typename Vector> 
  void AS_local_solve(using_qmr, const Matrix1 &A, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { qmr(A, x, b, P, iter); }

#ifdef GMM_USES_SUPERLU
  struct using_superlu {};

  template <typename P, typename Matrix>
  struct actual_precond<P, using_superlu, Matrix> {
    typedef typename linalg_traits<Matrix>::value_type value_type;
    typedef SuperLU_factor<value_type> APrecond;
    template <typename PR>
    static APrecond transform(const PR &) { return APrecond(); }
    static const APrecond &transform(const APrecond &PP) { return PP; }
  };

  template <typename Matrix1, typename Precond, typename Vector> 
  void AS_local_solve(using_superlu, const Matrix1 &, Vector &x,
		      const Vector &b, const Precond &P, iteration &iter)
  { P.solve(x, b); iter.set_iteration(1); }
#endif

  /* ******************************************************************** */
  /*		Additive Schwarz Linear system                            */
  /* ******************************************************************** */

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver>
  struct add_schwarz_mat{
    typedef typename linalg_traits<Matrix1>::value_type value_type;

    const Matrix1 *A;
    const std::vector<Matrix2> *vB;
    std::vector<Matrix2> vAloc;
    mutable iteration iter;
    double residu;
    mutable size_type itebilan;
    mutable std::vector<std::vector<value_type> > gi, fi;
    std::vector<typename actual_precond<Precond, local_solver,
					Matrix1>::APrecond> precond1;

    void init(const Matrix1 &A_, const std::vector<Matrix2> &vB_,
	      iteration iter_, const Precond &P, double residu_);

    add_schwarz_mat(void) {}
    add_schwarz_mat(const Matrix1 &A_, const std::vector<Matrix2> &vB_,
		iteration iter_, const Precond &P, double residu_)
    { init(A_, vB_, iter_, P, residu_); }
  };

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename local_solver>
  void add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver>::init(
       const Matrix1 &A_, const std::vector<Matrix2> &vB_,
       iteration iter_, const Precond &P, double residu_) {

    vB = &vB_; A = &A_; iter = iter_;
    residu = residu_;
    
    size_type nb_sub = vB->size();
    vAloc.resize(nb_sub);
    gi.resize(nb_sub); fi.resize(nb_sub);
    precond1.resize(nb_sub);
    std::fill(precond1.begin(), precond1.end(),
	      actual_precond<Precond, local_solver, Matrix1>::transform(P));
    itebilan = 0;
    
    if (iter.get_noisy()) cout << "Init pour sub dom ";
    for (size_type i = 0; i < nb_sub; ++i) {
      if (iter.get_noisy()) cout << i << " " << std::flush;
      Matrix2 Maux(mat_ncols((*vB)[i]), mat_nrows((*vB)[i]));
      
      gmm::resize(vAloc[i], mat_ncols((*vB)[i]), mat_ncols((*vB)[i]));      
      gmm::mult(gmm::transposed((*vB)[i]), *A, Maux);
      gmm::mult(Maux, (*vB)[i], vAloc[i]);
      precond1[i].build_with(vAloc[i]);
      gmm::resize(fi[i], mat_ncols((*vB)[i]));
      gmm::resize(gi[i], mat_ncols((*vB)[i]));
    }
    if (iter.get_noisy()) cout << "\n";
  }
  
  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename Vector2, typename Vector3, typename local_solver>
  void mult(const add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver> &M,
	    const Vector2 &p, Vector3 &q) {
    size_type itebilan = 0;
    mult(*(M.A), p, q);
    globaltolocal(q, M.fi, *(M.vB));
    for (size_type i = 0; i < M.fi.size(); ++i) {
      M.iter.init();
      AS_local_solve(local_solver(), (M.vAloc)[i], (M.gi)[i],
		     (M.fi)[i],(M.precond1)[i],M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }
    localtoglobal(M.gi, q, *(M.vB));
    if (M.iter.get_noisy() > 0) cout << "itebloc = " << itebilan << endl;
    M.itebilan += itebilan;
    M.iter.set_resmax((M.iter.get_resmax() + M.residu) * 0.5);
  }

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename Vector2, typename Vector3, typename local_solver>
  void mult(const add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver> &M,
	    const Vector2 &p, const Vector3 &q) {
    mult(M, p, const_cast<Vector3 &>(q));
  }

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename Vector2, typename Vector3, typename Vector4,
	    typename local_solver>
  void mult(const add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver> &M,
	    const Vector2 &p, const Vector3 &p2, Vector4 &q)
  { mult(M, p, q); add(p2, q); }

  template <typename Matrix1, typename Matrix2, typename Precond,
	    typename Vector2, typename Vector3, typename Vector4,
	    typename local_solver>
  void mult(const add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver> &M,
	    const Vector2 &p, const Vector3 &p2, const Vector4 &q)
  { mult(M, p, const_cast<Vector4 &>(q)); add(p2, q); }

  

  template <typename Matrix2, typename Vector2, typename Vector3>
  void globaltolocal(const Vector2 &f, std::vector<Vector3> &fi,
		       const std::vector<Matrix2> &vB) {
    for (size_type i = 0; i < fi.size(); ++i)
      gmm::mult(gmm::transposed(vB[i]), f, fi[i]);
  }

  template <typename Matrix2, typename Vector2, typename Vector3>
  void localtoglobal(const std::vector<Vector3> &fi, Vector2 &f, 
		     const std::vector<Matrix2> &vB) {
    gmm::clear(f);
    for (size_type i = 0; i < fi.size(); ++i) gmm::mult(vB[i], fi[i], f, f);
  }

  /* ******************************************************************** */
  /*		Additive Schwarz interfaced global solvers                */
  /* ******************************************************************** */

  template <typename ASM_type, typename Vect>
  void AS_global_solve(using_cg, const ASM_type &ASM, Vect &x,
		       const Vect &b, iteration &iter)
  { cg(ASM, x, b, *(ASM.A),  identity_matrix(), iter); }

  template <typename ASM_type, typename Vect>
  void AS_global_solve(using_gmres, const ASM_type &ASM, Vect &x,
		       const Vect &b, iteration &iter)
  { gmres(ASM, x, b, identity_matrix(), 100, iter); }

  template <typename ASM_type, typename Vect>
  void AS_global_solve(using_bicgstab, const ASM_type &ASM, Vect &x,
		       const Vect &b, iteration &iter)
  { bicgstab(ASM, x, b, identity_matrix(), iter); }

  template <typename ASM_type, typename Vect>
  void AS_global_solve(using_qmr,const ASM_type &ASM, Vect &x,
		       const Vect &b, iteration &iter)
  { qmr(ASM, x, b, identity_matrix(), iter); }

#ifdef GMM_USES_SUPERLU
  template <typename ASM_type, typename Vect>
  void AS_global_solve(using_superlu, const ASM_type &, Vect &,
		       const Vect &, iteration &) {
    DAL_THROW(failure_error,
     "You cannot use SuperLU as global solver in additive Schwarz meethod");
  }
#endif
  
  /* ******************************************************************** */
  /*		Sequential Linear Additive Schwarz method                 */
  /* ******************************************************************** */
  /* ref : Domain decomposition algorithms for the p-version finite       */
  /*       element method for elliptic problems, Luca F. Pavarino,        */
  /*       PhD thesis, Courant Institute of Mathematical Sciences, 1992.  */
  /* ******************************************************************** */

  /** Function to call if the ASM matrix is precomputed for successive solve
   * with the same system.
   */
  template <typename Matrix1, typename Matrix2,
	    typename Vector2, typename Vector3, typename Precond,
	    typename local_solver, typename global_solver>
  void sequential_additive_schwarz(
    add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver> &ASM, Vector3 &u,
    const Vector2 &f, iteration &iter, const global_solver&) {

    typedef typename linalg_traits<Matrix1>::value_type value_type;

    size_type nb_sub = ASM.vB->size(), nb_dof = gmm::vect_size(f);
    ASM.itebilan = 0;
    std::vector<value_type> g(nb_dof);

    for (size_type i = 0; i < nb_sub; ++i) {
      gmm::mult(gmm::transposed((*(ASM.vB))[i]), f, ASM.fi[i]);
      ASM.iter.init();
      AS_local_solve(local_solver(), ASM.vAloc[i], ASM.gi[i], ASM.fi[i],
		     ASM.precond1[i], ASM.iter);
      ASM.itebilan = std::max(ASM.itebilan, ASM.iter.get_iteration());
      gmm::mult((*(ASM.vB))[i], ASM.gi[i], g, g);
    }

    AS_global_solve(global_solver(), ASM, u, g, iter);
    if (iter.get_noisy())
      cout << "Total number of internal iterations : " << ASM.itebilan << endl;
  }

  /** Global function. Compute the ASM matrix and call the previous function.
      The ASM matrix represent the preconditionned linear system.
   */
  template <typename Matrix1, typename Matrix2,
	    typename Vector2, typename Vector3, typename Precond,
	    typename local_solver, typename global_solver>
  void sequential_additive_schwarz(const Matrix1 &A, Vector3 &u,
				  const Vector2 &f, const Precond &P,
				  const std::vector<Matrix2> &vB,
				  iteration &iter, local_solver,
				  global_solver) {
    iter.set_rhsnorm(vect_norm2(f));
    if (iter.get_rhsnorm() == 0.0) { gmm::clear(u); return; }
    iteration iter2 = iter; iter2.reduce_noisy();
    iter2.set_maxiter(size_type(-1));
    add_schwarz_mat<Matrix1, Matrix2, Precond, local_solver>
      ASM(A, vB, iter2, P, iter.get_resmax());
    sequential_additive_schwarz(ASM, u, f, iter, global_solver());
  }

  /* ******************************************************************** */
  /*		Sequential Non-Linear Additive Schwarz method             */
  /* ******************************************************************** */
  /* ref : Nonlinearly Preconditionned Inexact Newton Algorithms,         */
  /*       Xiao-Chuan Cai, David E. Keyes,                                */
  /*       SIAM J. Sci. Comp. 24: p183-200. .                             */
  /* ******************************************************************** */

  template <typename Matrixt, typename MatrixBi> 
  class NewtonAS_struct {
    
  public :
    typedef Matrixt tangent_matrix_type;
    typedef MatrixBi B_matrix_type;
    typedef typename linalg_traits<Matrixt>::value_type value_type;
    typedef std::vector<value_type> Vector;
    
    virtual size_type size(void) = 0;
    virtual const std::vector<MatrixBi> &get_vB() = 0;
    
    virtual void compute_F(Vector &f, Vector &x) = 0;
    virtual void compute_tangent_matrix(Matrixt &M, Vector &x) = 0;
    // compute Bi^T grad(F(X)) Bi
    virtual void compute_sub_tangent_matrix(Matrixt &Mloc, Vector &x,
					    size_type i) = 0;
    // compute Bi^T F(X)
    virtual void compute_sub_F(Vector &fi, Vector &x, size_type i) = 0;

    virtual ~NewtonAS_struct() {}
  };
  
  template <typename Matrixt, typename MatrixBi, typename Vector,
	    typename Precond, typename local_solver, typename global_solver>
  void Newton_additive_Schwarz(NewtonAS_struct<Matrixt, MatrixBi> &NS,
			       const Vector &u_,
			       iteration &iter, const Precond &P,
			       local_solver, global_solver) {
    Vector &u = const_cast<Vector &>(u_);
    typedef typename linalg_traits<Vector>::value_type value_type;
    typedef typename number_traits<value_type>::magnitude_type mtype;
    typedef actual_precond<Precond, local_solver, Matrixt> chgt_precond;
    
    double residu = iter.get_resmax();
    typename chgt_precond::APrecond PP = chgt_precond::transform(P);
    iter.set_rhsnorm(mtype(1));
    iteration iternc(iter);
    iternc.reduce_noisy(); iternc.set_maxiter(size_type(-1));
    iteration iter2(iternc);
    iteration iter3(iter2); iter3.reduce_noisy();
    iteration iter4(iter3);
    iternc.set_name("Local Newton");
    iter2.set_name("Linear System for Global Newton");
    iternc.set_resmax(residu/100.0);
    iter3.set_resmax(residu/1000.0);
    iter2.set_resmax(residu/100.0);
    iter4.set_resmax(residu/1000.0);
    std::vector<value_type> rhs(NS.size()), x(NS.size()), d(NS.size());
    std::vector<value_type> xi, xii, fi, di;
    Matrixt Mloc, M(NS.size(), NS.size());
    NS.compute_F(rhs, u);
    mtype act_res=gmm::vect_norm2(rhs), act_res_new(0), precond_res = act_res;
    mtype alpha, alpha_min=mtype(1)/mtype(16), alpha_mult=mtype(3)/mtype(4);
    mtype alpha_max_ratio(2);
    
    while(!iter.finished(std::min(act_res, precond_res))) {
      for (int SOR_step = 0;  SOR_step >= 0; --SOR_step) {
	gmm::clear(rhs);
	for (size_type isd = 0; isd < NS.get_vB().size(); ++isd) {
	  const MatrixBi &Bi = (NS.get_vB())[isd];
	  size_type si = mat_ncols(Bi);
	  gmm::resize(Mloc, si, si);
	  xi.resize(si); xii.resize(si); fi.resize(si); di.resize(si);
	  
	  iternc.init();
	  iternc.set_maxiter(30); // ?
	  if (iternc.get_noisy())
	    cout << "Non-linear local problem " << isd << endl;
	  gmm::clear(xi);
	  gmm::copy(u, x);
	  NS.compute_sub_F(fi, x, isd); gmm::scale(fi, value_type(-1));
	  mtype r = gmm::vect_norm2(fi), r_t;
	  if (r > value_type(0)) {
	    iternc.set_rhsnorm(std::max(r, mtype(1)));
	    while(!iternc.finished(r)) {
	      NS.compute_sub_tangent_matrix(Mloc, x, isd);
	      PP.build_with(Mloc);
	      iter3.init();
	      AS_local_solve(local_solver(), Mloc, di, fi, PP, iter3);
	      
	      for (alpha = mtype(1); alpha >= alpha_min; alpha *= alpha_mult) {
		gmm::add(xi, gmm::scaled(di, alpha), xii);
		gmm::mult(Bi, xii, u, x);
		NS.compute_sub_F(fi, x, isd); gmm::scale(fi, value_type(-1));
		if ((r_t = gmm::vect_norm2(fi)) <= r * alpha_max_ratio) break;
	      }
	      if (iternc.get_noisy()) cout << "(step=" << alpha << ")\t";
	      ++iternc; r = r_t; gmm::copy(xii, xi); 
	    }
	    if (SOR_step) gmm::mult(Bi, gmm::scaled(xii, 1.7), u, u);
	    gmm::mult(Bi, xii, rhs, rhs);
	  }
	}
	precond_res = gmm::vect_norm2(rhs);
	if (SOR_step) cout << "SOR step residu = " << precond_res << endl;
	if (precond_res < residu) break;
      }

      iter2.init();
      // solving linear system for the global Newton method
      NS.compute_tangent_matrix(M, u);
      add_schwarz_mat<Matrixt, MatrixBi, Precond, local_solver>
	ASM(M, NS.get_vB(), iter4, P, iter.get_resmax());
      AS_global_solve(global_solver(), ASM, d, rhs, iter2);

      for (alpha = mtype(1); alpha >= alpha_min; alpha *= alpha_mult) {
	gmm::add(gmm::scaled(d, alpha), u, x);
	NS.compute_F(rhs, x);
	act_res_new = gmm::vect_norm2(rhs);
	if (act_res_new <= act_res * alpha_max_ratio) break;
      }
      if (iter.get_noisy() > 1) cout << endl;
      act_res = act_res_new; 
      if (iter.get_noisy()) cout << "(step=" << alpha << ")\t unprecond res = " << act_res << " ";
      
      
      ++iter; gmm::copy(x, u);
    }
  }















  /* ******************************************************************** */
  /*		Old version, obsolete.                                    */
  /* ******************************************************************** */

  #define PRECOND ildltt_precond

  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI>
  struct schwarz_additif_matrix {
    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename dense_vector_type<value_type>::vector_type vector_type; 
    const Matrix1 *A;
    const std::vector<Matrix2> *ml1;
    const std::vector<Matrix3> *ml2;
    const std::vector<SUBI> *cor;
    mutable iteration iter;
    double residu;
    mutable size_t itebilan;
    std::vector<vector_type> *gi;
    std::vector<vector_type> *fi;
    
    std::vector<PRECOND<Matrix2> > *precond1;
    std::vector<PRECOND<Matrix3> > *precond2;

  };

  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3>
  int schwarz_additif(const Matrix1 &A, Vector3 &u,
		      const std::vector<Matrix2> &ml1,
		      const std::vector<Matrix3> &ml2,
		      const std::vector<SUBI> &cor,
		      const Vector2 &f, iteration &iter) IS_DEPRECATED;
  
  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3>
  int schwarz_additif(const Matrix1 &A, Vector3 &u,
		      const std::vector<Matrix2> &ml1,
		      const std::vector<Matrix3> &ml2,
		      const std::vector<SUBI> &cor,
		      const Vector2 &f, iteration &iter) {

    typedef typename linalg_traits<Matrix2>::value_type value_type;
    typedef typename dense_vector_type<value_type>::vector_type vector_type;

    size_type nb_sub = ml1.size() + ml2.size();
    size_t itebilan = 0;
    std::vector<vector_type> gi(nb_sub);
    std::vector<vector_type> fi(nb_sub);
    
    std::vector<PRECOND<Matrix2> > precond1(ml1.size());
    std::vector<PRECOND<Matrix3> > precond2(ml2.size());

    for (size_type i = 0; i < ml1.size(); ++i)
      precond1[i] = PRECOND<Matrix2>(ml1[i], 10, 1E-7);
    for (size_type i = 0; i < ml2.size(); ++i)
      precond2[i] = PRECOND<Matrix3>(ml2[i], 10, 1E-7);

    iter.set_rhsnorm(vect_norm2(f));
    if (iter.get_rhsnorm() == 0.0) { clear(u); return 0; }

    size_type ms = ml1.size();

    for (size_type i = 0; i < nb_sub; ++i) {
      size_type k = i < ms ? mat_nrows(ml1[i]) : mat_nrows(ml2[i-ms]);
      fi[i] = gi[i] = vector_type(k);
      clear(gi[i]);
    }

    size_type nb_dof = f.size();
    global_to_local(f, fi, cor);

    iteration iter2 = iter;
    iter2.reduce_noisy();

    for (size_type i = 0; i < ms; ++i) {
      iter2.init();
      cg(ml1[i], gi[i], fi[i], identity_matrix(), precond1[i], iter2);
      itebilan = std::max(itebilan, iter2.get_iteration());
    }
    for (size_type i = 0; i < ml2.size(); ++i) {
      iter2.init();
      cg(ml2[i], gi[i+ms], fi[i+ms],
	 identity_matrix(), precond2[i], iter2);
      itebilan = std::max(itebilan, iter2.get_iteration());
    }

    vector_type g(nb_dof);
    local_to_global(gi, g, cor);
    
    schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> SAM;
    SAM.A = &A; SAM.ml1 = &ml1; SAM.ml2 = &ml2; SAM.cor = &cor;
    SAM.precond1 = &precond1; SAM.precond2 = &precond2; 
    iter2.init();
    SAM.iter = iter2;
    SAM.residu = iter.get_resmax();
    // SAM.residu_act = 1E-2;
    SAM.gi = &gi; SAM.fi = &fi; SAM.itebilan = itebilan;
   
    cg(SAM, u, g, A, identity_matrix(), iter);

    return SAM.itebilan;
  }
  
  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, Vector3 &q) {

    size_type itebilan = 0;
    size_type ms = (M.ml1)->size();
    mult(*(M.A), p, q);
    global_to_local(q, *(M.fi), *(M.cor));
    for (size_type i = 0; i < (M.ml1)->size(); ++i) {
      M.iter.init();
      cg((*(M.ml1))[i], (*(M.gi))[i], (*(M.fi))[i],(*(M.precond1))[i], M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }

    for (size_type i = 0; i < (M.ml2)->size(); ++i) {
      M.iter.init();
      cg((*(M.ml2))[i],(*(M.gi))[i+ms], (*(M.fi))[i+ms], (*(M.precond2))[i],
	 M.iter);
      itebilan = std::max(itebilan, M.iter.get_iteration());
    }

    local_to_global(*(M.gi), q, *(M.cor));
    if (M.iter.get_noisy()) cout << "itebloc = " << itebilan << endl;
    M.itebilan += itebilan;
    M.iter.set_resmax((M.iter.get_resmax() + M.residu) * 0.5);
  }

  template <typename Matrix1, typename Matrix2, typename Matrix3,
	    typename SUBI,    typename Vector2, typename Vector3,
	    typename Vector4>
  void mult(const schwarz_additif_matrix<Matrix1, Matrix2, Matrix3, SUBI> &M,
	    const Vector2 &p, const Vector3 &p2, Vector4 &q)
  { mult(M, p, q); add(p2, q); }

  template <typename SUBI, typename Vector2, typename Vector3>
  void global_to_local(const Vector2 &f, std::vector<Vector3> &fi,
		       const std::vector<SUBI> &cor) {
    for (size_type i = 0; i < fi.size(); ++i) {
      typename linalg_traits<Vector3>::iterator it2 = fi[i].begin();
      for (size_type j = 0, l = cor[i].size(); j < l; ++j , ++it2)
        *it2 = f[cor[i].index(j)];
    }
  }

  template <typename SUBI, typename Vector2, typename Vector3>
  void local_to_global(const std::vector<Vector3> &fi, Vector2 &f, 
		       const std::vector<SUBI> &cor) {
    clear(f);
    for (size_type i = 0; i < fi.size(); ++i) {
      typename linalg_traits<Vector3>::const_iterator it2=fi[i].begin();
      for (size_type j = 0, l = cor[i].size(); j < l; ++j, ++it2) {
	f[cor[i].index(j)] += *it2;
      }
    }
  }


}


#endif //  GMM_SOLVERS_SCHWARZ_ADDITIVE_H__
