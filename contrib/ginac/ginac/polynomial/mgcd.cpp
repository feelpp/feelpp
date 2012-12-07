/** @file mgcd.cpp
 *
 *  Chinese remainder algorithm. */

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

#include "operators.h"
#include "chinrem_gcd.h"
#include "pgcd.h"
#include "collect_vargs.h"
#include "primes_factory.h"
#include "divide_in_z_p.h"
#include "poly_cra.h"
#include <numeric> // std::accumulate

#include <cln/integer.h>
#include <cln/integer_ring.h>
#include <cln/rational.h>
#include <cln/rational_ring.h>

namespace GiNaC {

static cln::cl_I extract_integer_content(ex& Apr, const ex& A)
{
	static const cln::cl_I n1(1);
	const numeric icont_ = A.integer_content();
	GINAC_ASSERT(cln::instanceof(icont_.to_cl_N(), cln::cl_RA_ring));
	if (cln::instanceof(icont_.to_cl_N(), cln::cl_I_ring)) {
		const cln::cl_I icont = cln::the<cln::cl_I>(icont_.to_cl_N());
		if (icont != 1) {
			Apr = (A/icont_).expand();
			return icont;
		} else {
			Apr = A;
			return n1;
		}
	} else {
		Apr = (A/icont_).expand();
		// A is a polynomail over rationals, so GCD is defined
		// up to arbitrary rational number.
		return n1;
	}
}

ex chinrem_gcd(const ex& A_, const ex& B_, const exvector& vars)
{
	ex A, B;
	const cln::cl_I a_icont = extract_integer_content(A, A_);
	const cln::cl_I b_icont = extract_integer_content(B, B_);
	const cln::cl_I c = cln::gcd(a_icont, b_icont);

	const cln::cl_I a_lc = integer_lcoeff(A, vars);
	const cln::cl_I b_lc = integer_lcoeff(B, vars);
	const cln::cl_I g_lc = cln::gcd(a_lc, b_lc);

	const ex& x(vars.back());
	exp_vector_t n = std::min(degree_vector(A, vars), degree_vector(B, vars));
	const int nTot = std::accumulate(n.begin(), n.end(), 0);
	const cln::cl_I A_max_coeff = to_cl_I(A.max_coefficient()); 
	const cln::cl_I B_max_coeff = to_cl_I(B.max_coefficient());

	const cln::cl_I lcoeff_limit = (cln::cl_I(1) << nTot)*cln::abs(g_lc)*
		std::min(A_max_coeff, B_max_coeff);


	cln::cl_I q = 0;
	ex H = 0;

	long p;
	primes_factory pfactory;
	while (true) {
		bool has_primes = pfactory(p, g_lc);
		if (!has_primes)
			throw chinrem_gcd_failed();

		const numeric pnum(p);
		ex Ap = A.smod(pnum);
		ex Bp = B.smod(pnum);
		ex Cp = pgcd(Ap, Bp, vars, p);

		const cln::cl_I g_lcp = smod(g_lc, p); 
		const cln::cl_I Cp_lc = integer_lcoeff(Cp, vars);
		const cln::cl_I nlc = smod(recip(Cp_lc, p)*g_lcp, p);
		Cp = (Cp*numeric(nlc)).expand().smod(pnum);
		exp_vector_t cp_deg = degree_vector(Cp, vars);
		if (zerop(cp_deg))
			return numeric(c);
		if (zerop(q)) {
			H = Cp;
			n = cp_deg;
			q = p;
		} else {
			if (cp_deg == n) {
				ex H_next = chinese_remainder(H, q, Cp, p);
				q = q*cln::cl_I(p);
				H = H_next;
			} else if (cp_deg < n) {
				// all previous homomorphisms are unlucky
				q = p;
				H = Cp;
				n = cp_deg;
			} else {
				// dp_deg > d_deg: current prime is bad
			}
		}
		if (q < lcoeff_limit)
			continue; // don't bother to do division checks
		ex C, dummy1, dummy2;
		extract_integer_content(C, H);
		if (divide_in_z_p(A, C, dummy1, vars, 0) && 
				divide_in_z_p(B, C, dummy2, vars, 0))
			return (numeric(c)*C).expand();
		// else: try more primes
	}
}

} // namespace GiNaC
