/** @file exam_clifford.cpp
 *
 *  Here we test GiNaC's Clifford algebra objects. */

/*
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "ginac.h"
using namespace GiNaC;

#include <iostream>
using namespace std;

const numeric half(1, 2);

static unsigned check_equal(const ex &e1, const ex &e2)
{
	ex e = normal(e1 - e2);
	if (!e.is_zero()) {
		clog << "(" << e1 << ") - (" << e2 << ") erroneously returned "
		     << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned check_equal_simplify(const ex &e1, const ex &e2)
{
	ex e = normal(simplify_indexed(e1) - e2);
	if (!e.is_zero()) {
		clog << "simplify_indexed(" << e1 << ") - (" << e2 << ") erroneously returned "
			 << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned check_equal_lst(const ex & e1, const ex & e2)
{
	for (unsigned int i = 0; i < e1.nops(); i++) {
		ex e = e1.op(i) - e2.op(i);
		if (!e.normal().is_zero()) {
			clog << "(" << e1 << ") - (" << e2 << ") erroneously returned "
			     << e << " instead of 0 (in the entry " << i  << ")" << endl;
			return 1;
		}
	}
	return 0;
}

static unsigned check_equal_simplify_term(const ex & e1, const ex & e2, idx & mu)
{
	ex e = expand_dummy_sum(normal(simplify_indexed(e1) - e2), true);

 	for (int j=0; j<4; j++) {
		ex esub = e.subs(
				is_a<varidx>(mu)
					? lst (
							mu == idx(j, mu.get_dim()),
							ex_to<varidx>(mu).toggle_variance() == idx(j, mu.get_dim())
						)
					: lst(mu == idx(j, mu.get_dim()))
			);
		if (!(canonicalize_clifford(esub).is_zero())) {
			clog << "simplify_indexed(" << e1 << ") - (" << e2 << ") erroneously returned "
				 << canonicalize_clifford(esub) << " instead of 0 for mu=" << j << endl;
			return 1;
		}
	}
	return 0;
}

static unsigned check_equal_simplify_term2(const ex & e1, const ex & e2)
{
 	ex e = expand_dummy_sum(normal(simplify_indexed(e1) - e2), true);
	if (!(canonicalize_clifford(e).is_zero())) {
		clog << "simplify_indexed(" << e1 << ") - (" << e2 << ") erroneously returned "
			 << canonicalize_clifford(e) << " instead of 0" << endl;
		return 1;
	}
	return 0;
}


static unsigned clifford_check1()
{
	// checks general identities and contractions

	unsigned result = 0;

	symbol dim("D");
	varidx mu(symbol("mu"), dim), nu(symbol("nu"), dim), rho(symbol("rho"), dim);
	ex e;

	e = dirac_ONE() * dirac_ONE();
	result += check_equal(e, dirac_ONE());

	e = dirac_ONE() * dirac_gamma(mu) * dirac_ONE();
	result += check_equal(e, dirac_gamma(mu));

	e = dirac_gamma(varidx(2, dim)) * dirac_gamma(varidx(1, dim)) *
	    dirac_gamma(varidx(1, dim)) * dirac_gamma(varidx(2, dim));
	result += check_equal(e, dirac_ONE());

	e = dirac_gamma(mu) * dirac_gamma(nu) *
	    dirac_gamma(nu.toggle_variance()) * dirac_gamma(mu.toggle_variance());
	result += check_equal_simplify(e, pow(dim, 2) * dirac_ONE());

	e = dirac_gamma(mu) * dirac_gamma(nu) *
	    dirac_gamma(mu.toggle_variance()) * dirac_gamma(nu.toggle_variance());
	result += check_equal_simplify(e, 2*dim*dirac_ONE()-pow(dim, 2)*dirac_ONE());

	e = dirac_gamma(nu.toggle_variance()) * dirac_gamma(rho.toggle_variance()) *
	    dirac_gamma(mu) * dirac_gamma(rho) * dirac_gamma(nu);
	e = e.simplify_indexed().collect(dirac_gamma(mu));
	result += check_equal(e, pow(2 - dim, 2).expand() * dirac_gamma(mu));

	return result;
}

static unsigned clifford_check2()
{
	// checks identities relating to gamma5

	unsigned result = 0;

	symbol dim("D");
	varidx mu(symbol("mu"), dim), nu(symbol("nu"), dim);
	ex e;

	e = dirac_gamma(mu) * dirac_gamma5() + dirac_gamma5() * dirac_gamma(mu);
	result += check_equal(e, 0);

	e = dirac_gamma5() * dirac_gamma(mu) * dirac_gamma5() + dirac_gamma(mu);
	result += check_equal(e, 0);

	return result;
}

static unsigned clifford_check3()
{
	// checks traces

	unsigned result = 0;

	symbol dim("D"), m("m"), q("q"), l("l"), ldotq("ldotq");
	varidx mu(symbol("mu"), dim), nu(symbol("nu"), dim), rho(symbol("rho"), dim),
	       sig(symbol("sig"), dim), kap(symbol("kap"), dim), lam(symbol("lam"), dim);
	ex e;

	e = dirac_gamma(mu);
	result += check_equal(dirac_trace(e), 0);

	e = dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho);
	result += check_equal(dirac_trace(e), 0);

	e = dirac_gamma5() * dirac_gamma(mu);
	result += check_equal(dirac_trace(e), 0);

	e = dirac_gamma5() * dirac_gamma(mu) * dirac_gamma(nu);
	result += check_equal(dirac_trace(e), 0);

	e = dirac_gamma5() * dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho);
	result += check_equal(dirac_trace(e), 0);

	scalar_products sp;
	sp.add(q, q, pow(q, 2));
	sp.add(l, l, pow(l, 2));
	sp.add(l, q, ldotq);

	e = pow(m, 2) * dirac_slash(q, dim) * dirac_slash(q, dim);
	e = dirac_trace(e).simplify_indexed(sp);
	result += check_equal(e, 4*pow(m, 2)*pow(q, 2));

	// cyclicity without gamma5
	e = dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig)
	  - dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig) * dirac_gamma(mu);
	e = dirac_trace(e);
	result += check_equal(e, 0);

	e = dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig) * dirac_gamma(kap) * dirac_gamma(lam)
	  - dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig) * dirac_gamma(kap) * dirac_gamma(lam) * dirac_gamma(mu);
	e = dirac_trace(e).expand();
	result += check_equal(e, 0);

	// cyclicity of gamma5 * S_4
	e = dirac_gamma5() * dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig)
	  - dirac_gamma(sig) * dirac_gamma5() * dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho);
	e = dirac_trace(e);
	result += check_equal(e, 0);

	// non-cyclicity of order D-4 of gamma5 * S_6
	e = dirac_gamma5() * dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig) * dirac_gamma(kap) * dirac_gamma(mu.toggle_variance())
	  + dim * dirac_gamma5() * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig) * dirac_gamma(kap);
	e = dirac_trace(e).simplify_indexed();
	e = (e / (dim - 4)).normal();
	result += check_equal(e, 8 * I * lorentz_eps(nu.replace_dim(4), rho.replace_dim(4), sig.replace_dim(4), kap.replace_dim(4)));

	// one-loop vacuum polarization in QED
	e = dirac_gamma(mu) *
	    (dirac_slash(l, dim) + dirac_slash(q, 4) + m * dirac_ONE()) *
	    dirac_gamma(mu.toggle_variance()) *
	    (dirac_slash(l, dim) + m * dirac_ONE());
	e = dirac_trace(e).simplify_indexed(sp);
	result += check_equal(e, 4*((2-dim)*l*l + (2-dim)*ldotq + dim*m*m).expand());

	e = dirac_slash(q, 4) *
	    (dirac_slash(l, dim) + dirac_slash(q, 4) + m * dirac_ONE()) *
	    dirac_slash(q, 4) *
	    (dirac_slash(l, dim) + m * dirac_ONE());
	e = dirac_trace(e).simplify_indexed(sp);
	result += check_equal(e, 4*(2*ldotq*ldotq + q*q*ldotq - q*q*l*l + q*q*m*m).expand());

	// stuff that had problems in the past
	ex prop = dirac_slash(q, dim) - m * dirac_ONE();
	e = dirac_slash(l, dim) * dirac_gamma5() * dirac_slash(l, dim) * prop;
	e = dirac_trace(dirac_slash(q, dim) * e) - dirac_trace(m * e)
	  - dirac_trace(prop * e);
	result += check_equal(e, 0);

	e = (dirac_gamma5() + dirac_ONE()) * dirac_gamma5();
	e = dirac_trace(e);
	result += check_equal(e, 4);

	// traces with multiple representation labels
	e = dirac_ONE(0) * dirac_ONE(1) / 16;
	result += check_equal(dirac_trace(e, 0), dirac_ONE(1) / 4);
	result += check_equal(dirac_trace(e, 1), dirac_ONE(0) / 4);
	result += check_equal(dirac_trace(e, 2), e);
	result += check_equal(dirac_trace(e, lst(0, 1)), 1);

	e = dirac_gamma(mu, 0) * dirac_gamma(mu.toggle_variance(), 1) * dirac_gamma(nu, 0) * dirac_gamma(nu.toggle_variance(), 1);
	result += check_equal_simplify(dirac_trace(e, 0), 4 * dim * dirac_ONE(1));
	result += check_equal_simplify(dirac_trace(e, 1), 4 * dim * dirac_ONE(0));
	// Fails with new tinfo mechanism because the order of gamme matrices with different rl depends on luck. 
	// TODO: better check.
	//result += check_equal_simplify(dirac_trace(e, 2), canonicalize_clifford(e)); // e will be canonicalized by the calculation of the trace
	result += check_equal_simplify(dirac_trace(e, lst(0, 1)), 16 * dim);

	return result;
}

static unsigned clifford_check4()
{
	// simplify_indexed()/dirac_trace() cross-checks

	unsigned result = 0;

	symbol dim("D");
	varidx mu(symbol("mu"), dim), nu(symbol("nu"), dim), rho(symbol("rho"), dim),
	       sig(symbol("sig"), dim), lam(symbol("lam"), dim);
	ex e, t1, t2;

	e = dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(mu.toggle_variance());
	t1 = dirac_trace(e).simplify_indexed();
	t2 = dirac_trace(e.simplify_indexed());
	result += check_equal((t1 - t2).expand(), 0);

	e = dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(sig) * dirac_gamma(mu.toggle_variance()) * dirac_gamma(lam);
	t1 = dirac_trace(e).simplify_indexed();
	t2 = dirac_trace(e.simplify_indexed());
	result += check_equal((t1 - t2).expand(), 0);

	e = dirac_gamma(sig) * dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(nu.toggle_variance()) * dirac_gamma(mu.toggle_variance());
	t1 = dirac_trace(e).simplify_indexed();
	t2 = dirac_trace(e.simplify_indexed());
	result += check_equal((t1 - t2).expand(), 0);

	e = dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(rho) * dirac_gamma(mu.toggle_variance()) * dirac_gamma(sig) * dirac_gamma(nu.toggle_variance());
	t1 = dirac_trace(e).simplify_indexed();
	t2 = dirac_trace(e.simplify_indexed());
	result += check_equal((t1 - t2).expand(), 0);

	return result;
}

static unsigned clifford_check5()
{
	// canonicalize_clifford() checks

	unsigned result = 0;

	symbol dim("D");
	varidx mu(symbol("mu"), dim), nu(symbol("nu"), dim), lam(symbol("lam"), dim);
	ex e;

	e = dirac_gamma(mu) * dirac_gamma(nu) + dirac_gamma(nu) * dirac_gamma(mu);
	result += check_equal(canonicalize_clifford(e), 2*dirac_ONE()*lorentz_g(mu, nu));

	e = (dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(lam)
	   + dirac_gamma(nu) * dirac_gamma(lam) * dirac_gamma(mu)
	   + dirac_gamma(lam) * dirac_gamma(mu) * dirac_gamma(nu)
	   - dirac_gamma(nu) * dirac_gamma(mu) * dirac_gamma(lam)
	   - dirac_gamma(lam) * dirac_gamma(nu) * dirac_gamma(mu)
	   - dirac_gamma(mu) * dirac_gamma(lam) * dirac_gamma(nu)) / 6
	  + lorentz_g(mu, nu) * dirac_gamma(lam)
	  - lorentz_g(mu, lam) * dirac_gamma(nu)
	  + lorentz_g(nu, lam) * dirac_gamma(mu)
	  - dirac_gamma(mu) * dirac_gamma(nu) * dirac_gamma(lam);
	result += check_equal(canonicalize_clifford(e), 0);

	return result;
}

/* We make two identical checks with metrics defined through a matrix in
 * the cases when used indexes have or have not variance.
 * To this end we recycle the code through the following macros */

template <typename IDX> unsigned clifford_check6(const matrix &A)
{
	unsigned result = 0;

	matrix A_symm(4,4), A2(4, 4);
	A_symm = A.add(A.transpose()).mul(half);
	A2 = A_symm.mul(A_symm);

	IDX v(symbol("v"), 4), nu(symbol("nu"), 4), mu(symbol("mu"), 4),
	       psi(symbol("psi"),4), lam(symbol("lambda"), 4),
	       xi(symbol("xi"), 4),  rho(symbol("rho"),4);
	ex mu_TOGGLE = is_a<varidx>(mu) ? ex_to<varidx>(mu).toggle_variance() : mu;
	ex nu_TOGGLE = is_a<varidx>(nu) ? ex_to<varidx>(nu).toggle_variance() : nu;
	ex rho_TOGGLE
		= is_a<varidx>(rho) ? ex_to<varidx>(rho).toggle_variance() : rho;

	ex e, e1;

/* checks general identities and contractions for clifford_unit*/
	e = dirac_ONE(2) * clifford_unit(mu, A, 2) * dirac_ONE(2);
	result += check_equal(e, clifford_unit(mu, A, 2));

	e = clifford_unit(IDX(2, 4), A) * clifford_unit(IDX(1, 4), A)
	  * clifford_unit(IDX(1, 4), A) * clifford_unit(IDX(2, 4), A);
	result += check_equal(e, A(1, 1) * A(2, 2) * dirac_ONE());

	e = clifford_unit(IDX(2, 4), A) * clifford_unit(IDX(1, 4), A)
	  * clifford_unit(IDX(1, 4), A) * clifford_unit(IDX(2, 4), A);
	result += check_equal(e, A(1, 1) * A(2, 2) * dirac_ONE());

	e = clifford_unit(nu, A) * clifford_unit(nu_TOGGLE, A);
	result += check_equal_simplify(e, A.trace() * dirac_ONE());

	e = clifford_unit(nu, A) * clifford_unit(nu, A);
	result += check_equal_simplify(e, indexed(A_symm, sy_symm(), nu, nu) * dirac_ONE());

	e = clifford_unit(nu, A) * clifford_unit(nu_TOGGLE, A) * clifford_unit(mu, A);
	result += check_equal_simplify(e, A.trace() * clifford_unit(mu, A));

	e = clifford_unit(nu, A) * clifford_unit(mu, A) * clifford_unit(nu_TOGGLE, A);
	
	result += check_equal_simplify_term(e,  2 * indexed(A_symm, sy_symm(), nu_TOGGLE, mu) *clifford_unit(nu, A)-A.trace()*clifford_unit(mu, A), mu);

	e = clifford_unit(nu, A) * clifford_unit(nu_TOGGLE, A)
	  * clifford_unit(mu, A) * clifford_unit(mu_TOGGLE, A);
	result += check_equal_simplify(e, pow(A.trace(), 2) * dirac_ONE());

	e = clifford_unit(mu, A) * clifford_unit(nu, A)
	  * clifford_unit(nu_TOGGLE, A) * clifford_unit(mu_TOGGLE, A);
	result += check_equal_simplify(e, pow(A.trace(), 2)  * dirac_ONE());

	e = clifford_unit(mu, A) * clifford_unit(nu, A)
	  * clifford_unit(mu_TOGGLE, A) * clifford_unit(nu_TOGGLE, A);

	result += check_equal_simplify_term2(e, 2*indexed(A_symm, sy_symm(), nu_TOGGLE, mu_TOGGLE) * clifford_unit(mu, A) * clifford_unit(nu, A) - pow(A.trace(), 2)*dirac_ONE());

	e = clifford_unit(mu_TOGGLE, A) * clifford_unit(nu, A)
	  * clifford_unit(mu, A) * clifford_unit(nu_TOGGLE, A);

	result += check_equal_simplify_term2(e, 2*indexed(A_symm, nu, mu) * clifford_unit(mu_TOGGLE, A) * clifford_unit(nu_TOGGLE, A) - pow(A.trace(), 2)*dirac_ONE());

	e = clifford_unit(nu_TOGGLE, A) * clifford_unit(rho_TOGGLE, A)
	  * clifford_unit(mu, A) * clifford_unit(rho, A) * clifford_unit(nu, A);
	e = e.simplify_indexed().collect(clifford_unit(mu, A));
	
	result += check_equal_simplify_term(e, 4* indexed(A_symm, sy_symm(), nu_TOGGLE,  rho)*indexed(A_symm, sy_symm(), rho_TOGGLE, mu) *clifford_unit(nu, A) 
	                                    - 2*A.trace() * (clifford_unit(rho, A) * indexed(A_symm, sy_symm(), rho_TOGGLE, mu) 
	                                                     + clifford_unit(nu, A) * indexed(A_symm, sy_symm(), nu_TOGGLE, mu)) + pow(A.trace(),2)* clifford_unit(mu, A), mu);

	e = clifford_unit(nu_TOGGLE, A) * clifford_unit(rho, A)
	  * clifford_unit(mu, A) * clifford_unit(rho_TOGGLE, A) * clifford_unit(nu, A);
	e = e.simplify_indexed().collect(clifford_unit(mu, A));
	
	result += check_equal_simplify_term(e, 4* indexed(A_symm, sy_symm(), nu_TOGGLE,  rho)*indexed(A_symm, sy_symm(), rho_TOGGLE, mu) *clifford_unit(nu, A) 
	                                    - 2*A.trace() * (clifford_unit(rho, A) * indexed(A_symm, sy_symm(), rho_TOGGLE, mu) 
	                                                     + clifford_unit(nu, A) * indexed(A_symm, sy_symm(), nu_TOGGLE, mu)) + pow(A.trace(),2)* clifford_unit(mu, A), mu);

	e = clifford_unit(mu, A) * clifford_unit(nu, A) + clifford_unit(nu, A) * clifford_unit(mu, A);
	result += check_equal(canonicalize_clifford(e), 2*dirac_ONE()*indexed(A_symm, sy_symm(), mu, nu));

	e = (clifford_unit(mu, A) * clifford_unit(nu, A) * clifford_unit(lam, A)
		 + clifford_unit(nu, A) * clifford_unit(lam, A) * clifford_unit(mu, A)
		 + clifford_unit(lam, A) * clifford_unit(mu, A) * clifford_unit(nu, A)
		 - clifford_unit(nu, A) * clifford_unit(mu, A) * clifford_unit(lam, A)
		 - clifford_unit(lam, A) * clifford_unit(nu, A) * clifford_unit(mu, A)
		 - clifford_unit(mu, A) * clifford_unit(lam, A) * clifford_unit(nu, A)) / 6
		+ indexed(A_symm, sy_symm(), mu, nu) * clifford_unit(lam, A)
		- indexed(A_symm, sy_symm(), mu, lam) * clifford_unit(nu, A)
		+ indexed(A_symm, sy_symm(), nu, lam) * clifford_unit(mu, A)
		- clifford_unit(mu, A) * clifford_unit(nu, A) * clifford_unit(lam, A);
	result += check_equal(canonicalize_clifford(e), 0);

/* lst_to_clifford() and clifford_inverse()  check*/
	realsymbol s("s"), t("t"), x("x"), y("y"), z("z");

	ex c = clifford_unit(nu, A, 1);
	e = lst_to_clifford(lst(t, x, y, z), mu, A, 1) * lst_to_clifford(lst(1, 2, 3, 4), c);
	e1 = clifford_inverse(e);
	result += check_equal_simplify_term2((e*e1).simplify_indexed(), dirac_ONE(1));

/* lst_to_clifford() and clifford_to_lst()  check for vectors*/
	e = lst(t, x, y, z);
	result += check_equal_lst(clifford_to_lst(lst_to_clifford(e, c), c, false), e);
	result += check_equal_lst(clifford_to_lst(lst_to_clifford(e, c), c, true), e);

/* lst_to_clifford() and clifford_to_lst()  check for pseudovectors*/
	e = lst(s, t, x, y, z);
	result += check_equal_lst(clifford_to_lst(lst_to_clifford(e, c), c, false), e);
	result += check_equal_lst(clifford_to_lst(lst_to_clifford(e, c), c, true), e);

/* Moebius map (both forms) checks for symmetric metrics only */
	matrix M1(2, 2),  M2(2, 2);
	c = clifford_unit(nu, A);

	e = clifford_moebius_map(0, dirac_ONE(), 
							 dirac_ONE(), 0, lst(t, x, y, z), A); 
/* this is just the inversion*/
	M1 = 0, dirac_ONE(),
		dirac_ONE(), 0;
	e1 = clifford_moebius_map(M1, lst(t, x, y, z), A); 
/* the inversion again*/
	result += check_equal_lst(e, e1);

	e1 = clifford_to_lst(clifford_inverse(lst_to_clifford(lst(t, x, y, z), mu, A)), c);
	result += check_equal_lst(e, e1);

	e = clifford_moebius_map(dirac_ONE(), lst_to_clifford(lst(1, 2, 3, 4), nu, A), 
							 0, dirac_ONE(), lst(t, x, y, z), A); 
/*this is just a shift*/
	M2 = dirac_ONE(), lst_to_clifford(lst(1, 2, 3, 4), c),
		0, dirac_ONE();
	e1 = clifford_moebius_map(M2, lst(t, x, y, z), c); 
/* the same shift*/
	result += check_equal_lst(e, e1);

	result += check_equal(e, lst(t+1, x+2, y+3, z+4));

/* Check the group law for Moebius maps */
	e = clifford_moebius_map(M1, ex_to<lst>(e1), c);
/*composition of M1 and M2*/
	e1 = clifford_moebius_map(M1.mul(M2), lst(t, x, y, z), c);
/* the product M1*M2*/
	result += check_equal_lst(e, e1);
	return result;
}

static unsigned clifford_check7(const ex & G, const symbol & dim)
{
	// checks general identities and contractions

	unsigned result = 0;

	varidx mu(symbol("mu"), dim), nu(symbol("nu"), dim), rho(symbol("rho"), dim),
	       psi(symbol("psi"),dim), lam(symbol("lambda"), dim), xi(symbol("xi"), dim);

	ex e;
	clifford unit = ex_to<clifford>(clifford_unit(mu, G));
	ex scalar = unit.get_metric(varidx(0, dim), varidx(0, dim));
	
	e = dirac_ONE() * dirac_ONE();
	result += check_equal(e, dirac_ONE());

	e = dirac_ONE() * clifford_unit(mu, G) * dirac_ONE();
	result += check_equal(e, clifford_unit(mu, G));

	e = clifford_unit(varidx(2, dim), G) * clifford_unit(varidx(1, dim), G)
	  * clifford_unit(varidx(1, dim), G) * clifford_unit(varidx(2, dim), G);
	result += check_equal(e, dirac_ONE()*pow(scalar, 2));

	e = clifford_unit(mu, G) * clifford_unit(nu, G)
	  * clifford_unit(nu.toggle_variance(), G) * clifford_unit(mu.toggle_variance(), G);
	result += check_equal_simplify(e, pow(dim*scalar, 2) * dirac_ONE());

	e = clifford_unit(mu, G) * clifford_unit(nu, G)
	  * clifford_unit(mu.toggle_variance(), G) * clifford_unit(nu.toggle_variance(), G);
	result += check_equal_simplify(e, (2*dim - pow(dim, 2))*pow(scalar,2)*dirac_ONE());

	e = clifford_unit(nu.toggle_variance(), G) * clifford_unit(rho.toggle_variance(), G)
	  * clifford_unit(mu, G) * clifford_unit(rho, G) * clifford_unit(nu, G);
	e = e.simplify_indexed().collect(clifford_unit(mu, G));
	result += check_equal(e, pow(scalar*(dim-2), 2).expand() * clifford_unit(mu, G));

	// canonicalize_clifford() checks, only for symmetric metrics
	if (is_a<indexed>(ex_to<clifford>(clifford_unit(mu, G)).get_metric()) &&
	    ex_to<symmetry>(ex_to<indexed>(ex_to<clifford>(clifford_unit(mu, G)).get_metric()).get_symmetry()).has_symmetry()) {
		e = clifford_unit(mu, G) * clifford_unit(nu, G) + clifford_unit(nu, G) * clifford_unit(mu, G);
		result += check_equal(canonicalize_clifford(e), 2*dirac_ONE()*unit.get_metric(nu, mu));
		
		e = (clifford_unit(mu, G) * clifford_unit(nu, G) * clifford_unit(lam, G)
			 + clifford_unit(nu, G) * clifford_unit(lam, G) * clifford_unit(mu, G)
			 + clifford_unit(lam, G) * clifford_unit(mu, G) * clifford_unit(nu, G)
			 - clifford_unit(nu, G) * clifford_unit(mu, G) * clifford_unit(lam, G)
			 - clifford_unit(lam, G) * clifford_unit(nu, G) * clifford_unit(mu, G)
			 - clifford_unit(mu, G) * clifford_unit(lam, G) * clifford_unit(nu, G)) / 6
			+ unit.get_metric(mu, nu) * clifford_unit(lam, G)
			- unit.get_metric(mu, lam) * clifford_unit(nu, G)
			+ unit.get_metric(nu, lam) * clifford_unit(mu, G)
			- clifford_unit(mu, G) * clifford_unit(nu, G) * clifford_unit(lam, G);
		result += check_equal(canonicalize_clifford(e), 0);
	} else {
		e = clifford_unit(mu, G) * clifford_unit(nu, G) + clifford_unit(nu, G) * clifford_unit(mu, G);
		result += check_equal(canonicalize_clifford(e), dirac_ONE()*(unit.get_metric(mu, nu) + unit.get_metric(nu, mu)));
		
		e = (clifford_unit(mu, G) * clifford_unit(nu, G) * clifford_unit(lam, G)
			 + clifford_unit(nu, G) * clifford_unit(lam, G) * clifford_unit(mu, G)
			 + clifford_unit(lam, G) * clifford_unit(mu, G) * clifford_unit(nu, G)
			 - clifford_unit(nu, G) * clifford_unit(mu, G) * clifford_unit(lam, G)
			 - clifford_unit(lam, G) * clifford_unit(nu, G) * clifford_unit(mu, G)
			 - clifford_unit(mu, G) * clifford_unit(lam, G) * clifford_unit(nu, G)) / 6
			+ half * (unit.get_metric(mu, nu) + unit.get_metric(nu, mu)) * clifford_unit(lam, G)
			- half * (unit.get_metric(mu, lam) + unit.get_metric(lam, mu)) * clifford_unit(nu, G)
			+ half * (unit.get_metric(nu, lam) + unit.get_metric(lam, nu)) * clifford_unit(mu, G)
			- clifford_unit(mu, G) * clifford_unit(nu, G) * clifford_unit(lam, G);
		result += check_equal(canonicalize_clifford(e), 0);
	}
	return result;
}

static unsigned clifford_check8()
{
	unsigned result = 0;

	realsymbol a("a");
	varidx mu(symbol("mu", "\\mu"), 1);

	ex e = clifford_unit(mu, diag_matrix(lst(-1))), e0 = e.subs(mu==0);
	result += ( exp(a*e0)*e0*e0 == -exp(e0*a) ) ? 0 : 1;

	return result;
}

unsigned exam_clifford()
{
	unsigned result = 0;
	
	cout << "examining clifford objects" << flush;

	result += clifford_check1(); cout << '.' << flush;
	result += clifford_check2(); cout << '.' << flush;
	result += clifford_check3(); cout << '.' << flush;
	result += clifford_check4(); cout << '.' << flush;
	result += clifford_check5(); cout << '.' << flush;

	// anticommuting, symmetric examples
	result += clifford_check6<varidx>(ex_to<matrix>(diag_matrix(lst(-1, 1, 1, 1))));
	result += clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-1, 1, 1, 1))));; cout << '.' << flush;
	result += clifford_check6<varidx>(ex_to<matrix>(diag_matrix(lst(-1, -1, -1, -1))))+clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-1, -1, -1, -1))));; cout << '.' << flush;
	result += clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-1, 1, 1, -1))))+clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-1, 1, 1, -1))));; cout << '.' << flush;
	result += clifford_check6<varidx>(ex_to<matrix>(diag_matrix(lst(-1, 0, 1, -1))))+clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-1, 0, 1, -1))));; cout << '.' << flush;
	result += clifford_check6<varidx>(ex_to<matrix>(diag_matrix(lst(-3, 0, 2, -1))))+clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-3, 0, 2, -1))));; cout << '.' << flush;

	realsymbol s("s"), t("t"); // symbolic entries in matric
	result += clifford_check6<varidx>(ex_to<matrix>(diag_matrix(lst(-1, 1, s, t))))+clifford_check6<idx>(ex_to<matrix>(diag_matrix(lst(-1, 1, s, t))));; cout << '.' << flush;

	matrix A(4, 4);
	A = 1, 0, 0, 0, // anticommuting, not symmetric, Tr=0
		0, -1, 0, 0,
		0, 0, 0, -1,
		0, 0, 1, 0; 
	result += clifford_check6<varidx>(A)+clifford_check6<idx>(A);; cout << '.' << flush;

	A = 1, 0, 0, 0, // anticommuting, not symmetric, Tr=2
		0, 1, 0, 0,
		0, 0, 0, -1,
		0, 0, 1, 0; 
	result += clifford_check6<varidx>(A)+clifford_check6<idx>(A);; cout << '.' << flush;

	A = 1, 0, 0, 0, // not anticommuting, symmetric, Tr=0
		0, -1, 0, 0,
		0, 0, 0, -1,
		0, 0, -1, 0; 
	result += clifford_check6<varidx>(A)+clifford_check6<idx>(A);; cout << '.' << flush;

	A = 1, 0, 0, 0, // not anticommuting, symmetric, Tr=2
		0, 1, 0, 0,
		0, 0, 0, -1,
		0, 0, -1, 0; 
	result += clifford_check6<varidx>(A)+clifford_check6<idx>(A);; cout << '.' << flush;

	A = 1, 1, 0, 0, // not anticommuting, not symmetric, Tr=4
		0, 1, 1, 0,
		0, 0, 1, 1,
		0, 0, 0, 1; 
	result += clifford_check6<varidx>(A)+clifford_check6<idx>(A);; cout << '.' << flush;

	symbol dim("D");
	result += clifford_check7(minkmetric(), dim); cout << '.' << flush;

	varidx chi(symbol("chi"), dim), xi(symbol("xi"), dim);
	result += clifford_check7(delta_tensor(xi, chi), dim); cout << '.' << flush;

	result += clifford_check7(lorentz_g(xi, chi), dim); cout << '.' << flush;

	result += clifford_check7(indexed(-2*minkmetric(), sy_symm(), xi, chi), dim); cout << '.' << flush;
	result += clifford_check7(-2*delta_tensor(xi, chi), dim); cout << '.' << flush;

	result += clifford_check8(); cout << '.' << flush;

	return result;
}

int main(int argc, char** argv)
{
	return exam_clifford();
}
