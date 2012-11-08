/** @file mod_gcd.cpp
 *
 *  Implementation of modular GCD. */

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

#include "upoly.h"
#include "gcd_euclid.h"
#include "cra_garner.h"
#include "debug.h"

#include <cln/numtheory.h>
#include <cln/random.h>
#include <stdexcept>

namespace GiNaC {

/**
 * @brief Remove the integer content from univariate polynomials A and B.
 *
 * As a byproduct compute the GCD of contents.
 */
static void remove_content(upoly& A, upoly& B, upoly::value_type& c)
{
	// remove the integer
	upoly::value_type acont, bcont;
	normalize_in_ring(A, &acont);
	normalize_in_ring(B, &bcont);
	c = gcd(acont, bcont);
};

/// Check if @a candidate divides both @a A and @a B
static bool
do_division_check(const upoly& A, const upoly& B, const upoly& candidate)
{
	upoly r1;
	remainder_in_ring(r1, A, candidate);
	if (r1.size() != 0)
		return false;

	upoly r2;
	remainder_in_ring(r2, B, candidate);
	if (r2.size() != 0)
		return false;

	return true;
}

/**
 * Given two GCD candidates H \in Z/q[x], and C \in Z/p[x] (where p is a prime)
 * compute the candidate in Z/(q*p) with CRA (chinise remainder algorithm)
 *
 * @param H \in Z/q[x] GCD candidate, will be updated by this function
 * @param q modulus of H, will NOT be updated by this function
 * @param C \in Z/p[x] GCD candidate, will be updated by this function
 * @param p modulus of C
 */
static void
update_the_candidate(upoly& H, const upoly::value_type& q,
	             const umodpoly& C,
	             const upoly::value_type& p,
		     const cln::cl_modint_ring& R)
{
	typedef upoly::value_type ring_t;
	std::vector<ring_t> moduli(2);
	moduli[0] = q;
	moduli[1] = p;
	if (H.size() < C.size())
		H.resize(C.size());

	for (std::size_t  i = C.size(); i-- != 0; ) {
		std::vector<ring_t> coeffs(2);
		coeffs[0] = H[i];
		coeffs[1] = R->retract(C[i]); 
		H[i] = integer_cra(coeffs, moduli);
	}
}

/// Convert Z/p[x] -> Z[x]
static void retract(upoly& p, const umodpoly& cp,
	            const cln::cl_modint_ring& Rp)
{
	p.resize(cp.size());
	for (std::size_t i = cp.size(); i-- != 0; )
		p[i] = Rp->retract(cp[i]);
}


/// Find the prime which is > p, and does NOT divide g
static void find_next_prime(cln::cl_I& p, const cln::cl_I& g)
{
	do {
		++p;
		p = nextprobprime(p);
	} while (zerop(mod(g, p)));
}

/// Compute the GCD of univariate polynomials A, B \in Z[x]
void mod_gcd(upoly& result, upoly A, upoly B)
{
	typedef upoly::value_type ring_t;
	ring_t content_gcd;
	// remove the integer content
	remove_content(A, B, content_gcd);

	// compute the coefficient bound for GCD(a, b)
	ring_t g = gcd(lcoeff(A), lcoeff(B));
	std::size_t max_gcd_degree = std::min(degree(A), degree(B));
	ring_t limit = (ring_t(1) << max_gcd_degree)*g*
		       std::min(max_coeff(A), max_coeff(B));
	ring_t q(0);
	upoly H;

	int count = 0;
	const ring_t p_threshold = ring_t(1) << (8*sizeof(void *));
	ring_t p = isqrt(std::min(max_coeff(A), max_coeff(B)));
	while (true) {
		if (count >= 8) {
			count = 0;
			if (p < p_threshold)
				p <<= 1;
			else
				p = p + (p >> 4);
		} else 
			++count;
		find_next_prime(p, g);

		// Map the polynomials onto Z/p[x]
		cln::cl_modint_ring Rp = cln::find_modint_ring(p);
		cln::cl_MI gp = Rp->canonhom(g);
		umodpoly ap(A.size()), bp(B.size());
		make_umodpoly(ap, A, Rp);
		make_umodpoly(bp, B, Rp);

		// Compute the GCD in Z/p[x]
		umodpoly cp;
		gcd_euclid(cp, ap, bp);
		bug_on(cp.size() == 0, "gcd(ap, bp) = 0, with ap = " <<
			                ap << ", and bp = " << bp);


		// Normalize the candidate so that its leading coefficient
		// is g mod p
		umodpoly::value_type norm_factor = gp*recip(lcoeff(cp));
		bug_on(zerop(norm_factor), "division in a field give 0");

		lcoeff(cp) = gp;
		for (std::size_t k = cp.size() - 1; k-- != 0; )
			cp[k] = cp[k]*norm_factor;


		// check for unlucky homomorphisms
		if (degree(cp) < max_gcd_degree) {
			q = p;
			max_gcd_degree = degree(cp);
			retract(H, cp, Rp);
		} else {
			update_the_candidate(H, q, cp, p, Rp);
			q = q*p;
		}

		if (q > limit) {
			result = H;
			normalize_in_ring(result);
			// if H divides both A and B it's a GCD
			if (do_division_check(A, B, result)) {
				result *= content_gcd;
				break; 
			}
			// H does not divide A and/or B, look for
			// another one
		} else if (degree(cp) == 0) {
			// Polynomials are relatively prime
			result.resize(1);
			result[0] = content_gcd;
			break;
		}
	}
}

} // namespace GiNaC
