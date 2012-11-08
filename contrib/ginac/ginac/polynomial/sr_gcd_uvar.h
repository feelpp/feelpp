/** @file sr_gcd_uvar.h
 *
 *  GCD function for univariate polynomials using PRS method. */

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

#ifndef GINAC_UPOLY_SR_GCD_H
#define GINAC_UPOLY_SR_GCD_H

#include "upoly.h"
#include "ring_traits.h"
#include "normalize.h"
#include "prem_uvar.h"

#include <limits>

namespace GiNaC {

/// Calculate GCD of two univariate polynomials @a a and @a b using
/// subresultant pseudo-remainder sequence method
template<typename T> static bool
sr_gcd_priv(T& g, T a, T b,
	    unsigned tries = std::numeric_limits<unsigned>::max())
{
	// handle zero polynomials
	if (a.empty()) {
		g.clear();
		return true;
	}
	if (b.empty()) {
		g.clear();
		return true;
	}

	typedef typename T::value_type ring_t;
	if (degree(a) < degree(b))
		swap(a, b);
	ring_t acont, bcont;
	normalize_in_ring(a, &acont);
	normalize_in_ring(b, &bcont);

	const ring_t one = get_ring_elt(b[0], 1);

	const ring_t gamma = gcd(acont, bcont);
	if (degree(b) == 0) {
		g.resize(1);
		g[0] = gamma;
		return true;
	}

	T r(std::min(degree(a), degree(b)));

	ring_t ri = one, psi = one;

	do {
		std::size_t delta = degree(a) - degree(b);
		pseudoremainder(r, a, b);
		if (r.empty()) {
			normalize_in_ring(b);
			b *= gamma;
			swap(b, g);
			return true;
		}
		a = b;

		ring_t ri_psi_delta = delta > 0 ?  ri*expt_pos(psi, delta) : ri;

		bool divisible_p = divide(b, r, ri_psi_delta);
		bug_on(!divisible_p, "division failed: r = " << r <<
			", ri = " << ri << ", psi = " << psi);

		if ((degree(b) == 0) && (degree(r) == 0)) {
			g.resize(1);
			g[0] = gamma;
			return true;
		}
		if (degree(b) == 0) {
			normalize_in_ring(r);
			r *= gamma;
			swap(r, g);
			return true;
		}

		ri = lcoeff(a);
		if (delta == 1)
			psi = ri;
		else if (delta) {
			const ring_t ri_delta = expt_pos(ri, delta);
			const ring_t psi_delta_1 = expt_pos(psi, delta - 1);
			bool sanity_check = div(psi, ri_delta, psi_delta_1);
			bug_on(!sanity_check, "division failed: ri = " << ri
				<< ", psi = " << psi << ", delta = " << delta);
		}
		if (--tries == 0)
			return false;
	} while (true);
}

} // namespace GiNaC

#endif // GINAC_UPOLY_SR_GCD_H
