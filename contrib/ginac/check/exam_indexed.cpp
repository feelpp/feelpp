/** @file exam_indexed.cpp
 *
 *  Here we test manipulations on GiNaC's indexed objects. */

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

static unsigned check_equal(const ex &e1, const ex &e2)
{
	ex e = e1 - e2;
	if (!e.is_zero()) {
		clog << e1 << "-" << e2 << " erroneously returned "
		     << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned check_equal_simplify(const ex &e1, const ex &e2)
{
	ex e = simplify_indexed(e1) - e2;
	if (!e.is_zero()) {
		clog << "simplify_indexed(" << e1 << ")-" << e2 << " erroneously returned "
		     << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned check_equal_simplify(const ex &e1, const ex &e2, const scalar_products &sp)
{
	ex e = simplify_indexed(e1, sp) - e2;
	if (!e.is_zero()) {
		clog << "simplify_indexed(" << e1 << ")-" << e2 << " erroneously returned "
		     << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned delta_check()
{
	// checks identities of the delta tensor

	unsigned result = 0;

	symbol s_i("i"), s_j("j"), s_k("k");
	idx i(s_i, 3), j(s_j, 3), k(s_k, 3);
	symbol A("A");

	// symmetry
	result += check_equal(delta_tensor(i, j), delta_tensor(j, i));

	// trace = dimension of index space
	result += check_equal(delta_tensor(i, i), 3);
	result += check_equal_simplify(delta_tensor(i, j) * delta_tensor(i, j), 3);

	// contraction with delta tensor
	result += check_equal_simplify(delta_tensor(i, j) * indexed(A, k), delta_tensor(i, j) * indexed(A, k));
	result += check_equal_simplify(delta_tensor(i, j) * indexed(A, j), indexed(A, i));
	result += check_equal_simplify(delta_tensor(i, j) * indexed(A, i), indexed(A, j));
	result += check_equal_simplify(delta_tensor(i, j) * delta_tensor(j, k) * indexed(A, i), indexed(A, k));

	return result;
}

static unsigned metric_check()
{
	// checks identities of the metric tensor

	unsigned result = 0;

	symbol s_mu("mu"), s_nu("nu"), s_rho("rho"), s_sigma("sigma");
	varidx mu(s_mu, 4), nu(s_nu, 4), rho(s_rho, 4), sigma(s_sigma, 4);
	symbol A("A");

	// becomes delta tensor if indices have opposite variance
	result += check_equal(metric_tensor(mu, nu.toggle_variance()), delta_tensor(mu, nu.toggle_variance()));

	// scalar contraction = dimension of index space
	result += check_equal(metric_tensor(mu, mu.toggle_variance()), 4);
	result += check_equal_simplify(metric_tensor(mu, nu) * metric_tensor(mu.toggle_variance(), nu.toggle_variance()), 4);

	// contraction with metric tensor
	result += check_equal_simplify(metric_tensor(mu, nu) * indexed(A, nu), metric_tensor(mu, nu) * indexed(A, nu));
	result += check_equal_simplify(metric_tensor(mu, nu) * indexed(A, nu.toggle_variance()), indexed(A, mu));
	result += check_equal_simplify(metric_tensor(mu, nu) * indexed(A, mu.toggle_variance()), indexed(A, nu));
	result += check_equal_simplify(metric_tensor(mu, nu) * metric_tensor(mu.toggle_variance(), rho.toggle_variance()) * indexed(A, nu.toggle_variance()), indexed(A, rho.toggle_variance()));
	result += check_equal_simplify(metric_tensor(mu, rho) * metric_tensor(nu, sigma) * indexed(A, rho.toggle_variance(), sigma.toggle_variance()), indexed(A, mu, nu));
	result += check_equal_simplify(indexed(A, mu.toggle_variance()) * metric_tensor(mu, nu) - indexed(A, mu.toggle_variance()) * metric_tensor(nu, mu), 0);
	result += check_equal_simplify(indexed(A, mu.toggle_variance(), nu.toggle_variance()) * metric_tensor(nu, rho), indexed(A, mu.toggle_variance(), rho));

	// contraction with delta tensor yields a metric tensor
	result += check_equal_simplify(delta_tensor(mu, nu.toggle_variance()) * metric_tensor(nu, rho), metric_tensor(mu, rho));
	result += check_equal_simplify(metric_tensor(mu, nu) * indexed(A, nu.toggle_variance()) * delta_tensor(mu.toggle_variance(), rho), indexed(A, rho));

	return result;
}

static unsigned epsilon_check()
{
	// checks identities of the epsilon tensor

	unsigned result = 0;

	symbol s_mu("mu"), s_nu("nu"), s_rho("rho"), s_sigma("sigma"), s_tau("tau");
	symbol d("d");
	varidx mu(s_mu, 4), nu(s_nu, 4), rho(s_rho, 4), sigma(s_sigma, 4), tau(s_tau, 4);
	varidx mu_co(s_mu, 4, true), nu_co(s_nu, 4, true), rho_co(s_rho, 4, true), sigma_co(s_sigma, 4, true), tau_co(s_tau, 4, true);

	// antisymmetry
	result += check_equal(lorentz_eps(mu, nu, rho, sigma) + lorentz_eps(sigma, rho, mu, nu), 0);

	// convolution is zero
	result += check_equal(lorentz_eps(mu, nu, rho, nu_co), 0);
	result += check_equal(lorentz_eps(mu, nu, mu_co, nu_co), 0);
	result += check_equal_simplify(lorentz_g(mu_co, nu_co) * lorentz_eps(mu, nu, rho, sigma), 0);

	// contraction with symmetric tensor is zero
	result += check_equal_simplify(lorentz_eps(mu, nu, rho, sigma) * indexed(d, sy_symm(), mu_co, nu_co), 0);
	result += check_equal_simplify(lorentz_eps(mu, nu, rho, sigma) * indexed(d, sy_symm(), nu_co, sigma_co, rho_co), 0);
	result += check_equal_simplify(lorentz_eps(mu, nu, rho, sigma) * indexed(d, mu_co) * indexed(d, nu_co), 0);
	result += check_equal_simplify(lorentz_eps(mu_co, nu, rho, sigma) * indexed(d, mu) * indexed(d, nu_co), 0);
	ex e = lorentz_eps(mu, nu, rho, sigma) * indexed(d, mu_co) - lorentz_eps(mu_co, nu, rho, sigma) * indexed(d, mu);
	result += check_equal_simplify(e, 0);

	// contractions of epsilon tensors
	result += check_equal_simplify(lorentz_eps(mu, nu, rho, sigma) * lorentz_eps(mu_co, nu_co, rho_co, sigma_co), -24);
	result += check_equal_simplify(lorentz_eps(tau, nu, rho, sigma) * lorentz_eps(mu_co, nu_co, rho_co, sigma_co), -6 * delta_tensor(tau, mu_co));

	return result;
}

DECLARE_FUNCTION_2P(symm_fcn)
REGISTER_FUNCTION(symm_fcn, set_symmetry(sy_symm(0, 1)));
DECLARE_FUNCTION_2P(anti_fcn)
REGISTER_FUNCTION(anti_fcn, set_symmetry(sy_anti(0, 1)));

static unsigned symmetry_check()
{
	// check symmetric/antisymmetric objects

	unsigned result = 0;

	idx i(symbol("i"), 3), j(symbol("j"), 3), k(symbol("k"), 3), l(symbol("l"), 3);
	symbol A("A"), B("B"), C("C");
	ex e;

	result += check_equal(indexed(A, sy_symm(), i, j), indexed(A, sy_symm(), j, i));
	result += check_equal(indexed(A, sy_anti(), i, j) + indexed(A, sy_anti(), j, i), 0);
	result += check_equal(indexed(A, sy_anti(), i, j, k) - indexed(A, sy_anti(), j, k, i), 0);
	e = indexed(A, sy_symm(), i, j, k) *
	    indexed(B, sy_anti(), l, k, i);
	result += check_equal_simplify(e, 0);
	e = indexed(A, sy_symm(), i, i, j, j) *
	    indexed(B, sy_anti(), k, l); // GiNaC 0.8.0 had a bug here
	result += check_equal_simplify(e, e);

	symmetry R = sy_symm(sy_anti(0, 1), sy_anti(2, 3));
	e = indexed(A, R, i, j, k, l) + indexed(A, R, j, i, k, l);
	result += check_equal(e, 0);
	e = indexed(A, R, i, j, k, l) + indexed(A, R, i, j, l, k);
	result += check_equal(e, 0);
	e = indexed(A, R, i, j, k, l) - indexed(A, R, j, i, l, k);
	result += check_equal(e, 0);
	e = indexed(A, R, i, j, k, l) + indexed(A, R, k, l, j, i);
	result += check_equal(e, 0);

	e = indexed(A, i, j);
	result += check_equal(symmetrize(e) + antisymmetrize(e), e);
	e = indexed(A, sy_symm(), i, j, k, l);
	result += check_equal(symmetrize(e), e);
	result += check_equal(antisymmetrize(e), 0);
	e = indexed(A, sy_anti(), i, j, k, l);
	result += check_equal(symmetrize(e), 0);
	result += check_equal(antisymmetrize(e), e);

	e = (indexed(A, sy_anti(), i, j, k, l) * (indexed(B, j) * indexed(C, k) + indexed(B, k) * indexed(C, j)) + indexed(B, i, l)).expand();
	result += check_equal_simplify(e, indexed(B, i, l));

	result += check_equal(symm_fcn(0, 1) + symm_fcn(1, 0), 2*symm_fcn(0, 1));
	result += check_equal(anti_fcn(0, 1) + anti_fcn(1, 0), 0);
	result += check_equal(anti_fcn(0, 0), 0);

	return result;
}

static unsigned scalar_product_check()
{
	// check scalar product replacement

	unsigned result = 0;

    idx i(symbol("i"), 3), j(symbol("j"), 3);
    symbol A("A"), B("B"), C("C");
	ex e;

    scalar_products sp;
    sp.add(A, B, 0); // A and B are orthogonal
    sp.add(A, C, 0); // A and C are orthogonal
    sp.add(A, A, 4); // A^2 = 4 (A has length 2)

    e = (indexed(A + B, i) * indexed(A + C, i)).expand(expand_options::expand_indexed);
	result += check_equal_simplify(e, indexed(B, i) * indexed(C, i) + 4, sp);
	e = indexed(A, i, i) * indexed(B, j, j); // GiNaC 0.8.0 had a bug here
	result += check_equal_simplify(e, e, sp);

	return result;
}

static unsigned edyn_check()
{
	// Relativistic electrodynamics

	// Test 1: check transformation laws of electric and magnetic fields by
	// applying a Lorentz boost to the field tensor

	unsigned result = 0;

	symbol beta("beta");
	ex gamma = 1 / sqrt(1 - pow(beta, 2));
	symbol Ex("Ex"), Ey("Ey"), Ez("Ez");
	symbol Bx("Bx"), By("By"), Bz("Bz");

	// Lorentz transformation matrix (boost along x axis)
	matrix L(4, 4);
	L =       gamma, -beta*gamma, 0, 0,
	    -beta*gamma,       gamma, 0, 0,
	              0,           0, 1, 0,
	              0,           0, 0, 1;

	// Electromagnetic field tensor
	matrix F(4, 4);
	F =  0, -Ex, -Ey, -Ez,
		Ex,   0, -Bz,  By,
		Ey,  Bz,   0, -Bx,
		Ez, -By,  Bx,   0;

	// Indices
	symbol s_mu("mu"), s_nu("nu"), s_rho("rho"), s_sigma("sigma");
	varidx mu(s_mu, 4), nu(s_nu, 4), rho(s_rho, 4), sigma(s_sigma, 4);

	// Apply transformation law of second rank tensor
	ex e = (indexed(L, mu, rho.toggle_variance())
	      * indexed(L, nu, sigma.toggle_variance())
	      * indexed(F, rho, sigma)).simplify_indexed();

	// Extract transformed electric and magnetic fields
	ex Ex_p = e.subs(lst(mu == 1, nu == 0)).normal();
	ex Ey_p = e.subs(lst(mu == 2, nu == 0)).normal();
	ex Ez_p = e.subs(lst(mu == 3, nu == 0)).normal();
	ex Bx_p = e.subs(lst(mu == 3, nu == 2)).normal();
	ex By_p = e.subs(lst(mu == 1, nu == 3)).normal();
	ex Bz_p = e.subs(lst(mu == 2, nu == 1)).normal();

	// Check results
	result += check_equal(Ex_p, Ex);
	result += check_equal(Ey_p, gamma * (Ey - beta * Bz));
	result += check_equal(Ez_p, gamma * (Ez + beta * By));
	result += check_equal(Bx_p, Bx);
	result += check_equal(By_p, gamma * (By + beta * Ez));
	result += check_equal(Bz_p, gamma * (Bz - beta * Ey));

	// Test 2: check energy density and Poynting vector of electromagnetic field

	// Minkowski metric
	ex eta = diag_matrix(lst(1, -1, -1, -1));

	// Covariant field tensor
	ex F_mu_nu = (indexed(eta, mu.toggle_variance(), rho.toggle_variance())
	            * indexed(eta, nu.toggle_variance(), sigma.toggle_variance())
	            * indexed(F, rho, sigma)).simplify_indexed();

	// Energy-momentum tensor
	ex T = (-indexed(eta, rho, sigma) * F_mu_nu.subs(s_nu == s_rho) 
	        * F_mu_nu.subs(lst(s_mu == s_nu, s_nu == s_sigma))
	      + indexed(eta, mu.toggle_variance(), nu.toggle_variance())
	        * F_mu_nu.subs(lst(s_mu == s_rho, s_nu == s_sigma))
	        * indexed(F, rho, sigma) / 4).simplify_indexed() / (4 * Pi);

	// Extract energy density and Poynting vector
	ex E = T.subs(lst(s_mu == 0, s_nu == 0)).normal();
	ex Px = T.subs(lst(s_mu == 0, s_nu == 1));
	ex Py = T.subs(lst(s_mu == 0, s_nu == 2)); 
	ex Pz = T.subs(lst(s_mu == 0, s_nu == 3));

	// Check results
	result += check_equal(E, (Ex*Ex+Ey*Ey+Ez*Ez+Bx*Bx+By*By+Bz*Bz) / (8 * Pi));
	result += check_equal(Px, (Ez*By-Ey*Bz) / (4 * Pi));
	result += check_equal(Py, (Ex*Bz-Ez*Bx) / (4 * Pi));
	result += check_equal(Pz, (Ey*Bx-Ex*By) / (4 * Pi));

	return result;
}

static unsigned spinor_check()
{
	// check identities of the spinor metric

	unsigned result = 0;

	symbol psi("psi");
	spinidx A(symbol("A")), B(symbol("B")), C(symbol("C")), D(symbol("D"));
	ex A_co = A.toggle_variance(), B_co = B.toggle_variance();
	ex e;

	e = spinor_metric(A_co, B_co) * spinor_metric(A, B);
	result += check_equal_simplify(e, 2);
	e = spinor_metric(A_co, B_co) * spinor_metric(B, A);
	result += check_equal_simplify(e, -2);
	e = spinor_metric(A_co, B_co) * spinor_metric(A, C);
	result += check_equal_simplify(e, delta_tensor(B_co, C));
	e = spinor_metric(A_co, B_co) * spinor_metric(B, C);
	result += check_equal_simplify(e, -delta_tensor(A_co, C));
	e = spinor_metric(A_co, B_co) * spinor_metric(C, A);
	result += check_equal_simplify(e, -delta_tensor(B_co, C));
	e = spinor_metric(A, B) * indexed(psi, B_co);
	result += check_equal_simplify(e, indexed(psi, A));
	e = spinor_metric(A, B) * indexed(psi, A_co);
	result += check_equal_simplify(e, -indexed(psi, B));
	e = spinor_metric(A_co, B_co) * indexed(psi, B);
	result += check_equal_simplify(e, -indexed(psi, A_co));
	e = spinor_metric(A_co, B_co) * indexed(psi, A);
	result += check_equal_simplify(e, indexed(psi, B_co));
	e = spinor_metric(D, A) * spinor_metric(A_co, B_co) * spinor_metric(B, C) - spinor_metric(D, A_co) * spinor_metric(A, B_co) * spinor_metric(B, C);
	result += check_equal_simplify(e, 0);

	return result;
}

static unsigned dummy_check()
{
	// check dummy index renaming/repositioning

	unsigned result = 0;

	symbol p("p"), q("q");
	idx i(symbol("i"), 3), j(symbol("j"), 3), n(symbol("n"), 3);
	varidx mu(symbol("mu"), 4), nu(symbol("nu"), 4);
	ex e;

	e = indexed(p, i) * indexed(q, i) - indexed(p, j) * indexed(q, j);
	result += check_equal_simplify(e, 0);

	e = indexed(p, i) * indexed(p, i) * indexed(q, j) * indexed(q, j)
	  - indexed(p, n) * indexed(p, n) * indexed(q, j) * indexed(q, j);
	result += check_equal_simplify(e, 0);

	e = indexed(p, mu, mu.toggle_variance()) - indexed(p, nu, nu.toggle_variance());
	result += check_equal_simplify(e, 0);

	e = indexed(p, mu.toggle_variance(), nu, mu) * indexed(q, i)
	  - indexed(p, mu, nu, mu.toggle_variance()) * indexed(q, i);
	result += check_equal_simplify(e, 0);

	e = indexed(p, mu, mu.toggle_variance()) - indexed(p, nu.toggle_variance(), nu);
	result += check_equal_simplify(e, 0);
	e = indexed(p, mu.toggle_variance(), mu) - indexed(p, nu, nu.toggle_variance());
	result += check_equal_simplify(e, 0);

	// GiNaC 1.2.1 had a bug here because p.i*p.i -> (p.i)^2
	e = indexed(p, i) * indexed(p, i) * indexed(p, j) + indexed(p, j);
	ex fi = exprseq(e.get_free_indices());
	if (!fi.is_equal(exprseq(j))) {
		clog << "get_free_indices(" << e << ") erroneously returned "
		     << fi << " instead of (.j)" << endl;
		++result;
	}

	return result;
}

unsigned exam_indexed()
{
	unsigned result = 0;
	
	cout << "examining indexed objects" << flush;

	result += delta_check();  cout << '.' << flush;
	result += metric_check();  cout << '.' << flush;
	result += epsilon_check();  cout << '.' << flush;
	result += symmetry_check();  cout << '.' << flush;
	result += scalar_product_check();  cout << '.' << flush;
	result += edyn_check();  cout << '.' << flush;
	result += spinor_check(); cout << '.' << flush;
	result += dummy_check(); cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_indexed();
}
