/** @file pgcd.cpp
 *
 *  GCD for polynomials in prime field. */

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

#include "pgcd.h"
#include "collect_vargs.h"
#include "smod_helpers.h"
#include "euclid_gcd_wrap.h"
#include "eval_point_finder.h"
#include "newton_interpolate.h"
#include "divide_in_z_p.h"

namespace GiNaC {

extern void
primpart_content(ex& pp, ex& c, ex e, const exvector& vars, const long p);

// Computes the GCD of two polynomials over a prime field.
// Based on Algorithm 7.2 from "Algorithms for Computer Algebra"
// A and B are considered as Z_p[x_n][x_0, \ldots, x_{n-1}], that is,
// as a polynomials in variables x_0, \ldots x_{n-1} having coefficients
// from the ring Z_p[x_n]
ex pgcd(const ex& A, const ex& B, const exvector& vars, const long p)
{
	static const ex ex1(1);
	if (A.is_zero())
		return B;

	if (B.is_zero())
		return A;

	if (is_a<numeric>(A))
		return ex1;

	if (is_a<numeric>(B))
		return ex1;

	// Checks for univariate polynomial
	if (vars.size() == 1) {
		ex ret = euclid_gcd(A, B, vars[0], p); // Univariate GCD
		return ret;
	}
	const ex& mainvar(vars.back());

	// gcd of the contents
	ex H = 0, Hprev = 0; // GCD candidate
	ex newton_poly = 1;  // for Newton Interpolation

	// Contents and primparts of A and B
	ex Aprim, contA;
	primpart_content(Aprim, contA, A, vars, p);
	ex Bprim, contB;
	primpart_content(Bprim, contB, B, vars, p);
	// gcd of univariate polynomials
	const ex cont_gcd = euclid_gcd(contA, contB, mainvar, p);

	exvector restvars = vars;
	restvars.pop_back();
	const ex AL = lcoeff_wrt(Aprim, restvars);
	const ex BL = lcoeff_wrt(Bprim, restvars);
	// gcd of univariate polynomials
	const ex lc_gcd = euclid_gcd(AL, BL, mainvar, p);

	// The estimate of degree of the gcd of Ab and Bb
	exp_vector_t gcd_deg = std::min(degree_vector(Aprim, restvars),
					degree_vector(Bprim, restvars));
	eval_point_finder::value_type b;

	eval_point_finder find_eval_point(p);
	const numeric pn(p);
	do {
		// Find a `good' evaluation point b.
		bool has_more_pts = find_eval_point(b, lc_gcd, mainvar);
		// If there are no more possible evaluation points, bail out
		if (!has_more_pts)
			throw pgcd_failed();

		const numeric bn(b);
		// Evaluate the polynomials in b
		ex Ab = Aprim.subs(mainvar == bn).smod(pn);
		ex Bb = Bprim.subs(mainvar == bn).smod(pn);
		ex Cb = pgcd(Ab, Bb, restvars, p);

		// Set the correct the leading coefficient
		const cln::cl_I lcb_gcd =
			smod(to_cl_I(lc_gcd.subs(mainvar == bn)), p);
		const cln::cl_I Cblc = integer_lcoeff(Cb, restvars);
		const cln::cl_I correct_lc = smod(lcb_gcd*recip(Cblc, p), p);
		Cb = (Cb*numeric(correct_lc)).smod(pn);

		const exp_vector_t img_gcd_deg = degree_vector(Cb, restvars);
		// Test for relatively prime polynomials
		if (zerop(img_gcd_deg))
			return cont_gcd;
		// Test for unlucky homomorphisms
		if (img_gcd_deg < gcd_deg) {
			// The degree decreased, previous homomorphisms were
			// bad, so we have to start it all over.
			H = Cb;
			newton_poly = mainvar - numeric(b);
			Hprev = 0;
			gcd_deg  = img_gcd_deg;
			continue;
		} 
		if (img_gcd_deg > gcd_deg) {
			// The degree of images GCD is too high, this
			// evaluation point is bad. Skip it.
			continue;
		}

		// Image has the same degree as the previous one
		// (or at least not higher than the limit)
		Hprev = H;
		H = newton_interp(Cb, b, H, newton_poly, mainvar, p);
		newton_poly = newton_poly*(mainvar - b);

		// try to reduce the number of division tests.
		const ex H_lcoeff = lcoeff_wrt(H, restvars);

		if (H_lcoeff.is_equal(lc_gcd)) {
			ex C /* primitive part of H */, contH /* dummy */;
			primpart_content(C, contH, H, vars, p);
			// Normalize GCD so that leading coefficient is 1
			const cln::cl_I Clc = recip(integer_lcoeff(C, vars), p);
			C = (C*numeric(Clc)).expand().smod(pn);

			ex dummy1, dummy2;

			if (divide_in_z_p(Aprim, C, dummy1, vars, p) &&
					divide_in_z_p(Bprim, C, dummy2, vars, p))
				return (cont_gcd*C).expand().smod(pn);
			// else continue building the candidate
		} 
	} while(true);
	throw pgcd_failed();
}

} // namespace GiNaC
