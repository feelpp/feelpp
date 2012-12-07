/** @file heur_gcd_uvar.h
 *
 *  Heuristic GCD code. */

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

#ifndef GINAC_UPOLY_HEUR_GCD_H
#define GINAC_UPOLY_HEUR_GCD_H

#include "upoly.h"
#include "ring_traits.h"
#include "normalize.h"
#include "remainder.h"
#include "eval_uvar.h"
#include "interpolate_padic_uvar.h"

#include <algorithm>

namespace GiNaC {

/// Compute GCD of primitive univariate polynomials.
template<typename T> static bool
heur_gcd_z_pp(T& g, const T& a, const T& b, unsigned max_tries = 66)
{
	typedef typename T::value_type ring_t;
	const ring_t n73794 = get_ring_elt(b[0], 73794);
	const ring_t n27011 = get_ring_elt(b[0], 27011);
	const std::size_t maxdeg = std::max(degree(a), degree(b));
	T r, gg;
	gg.reserve(maxdeg);
	r.reserve(maxdeg);

	// find the evaluation point
	ring_t xi = (std::min(max_coeff(a), max_coeff(b)) + 1) << 1;

	do {
		const ring_t av = eval(a, xi);
		const ring_t bv = eval(b, xi); 
		const ring_t gamma = gcd(av, bv);
		interpolate(gg, gamma, xi, maxdeg);
		normalize_in_ring(gg);
		remainder_in_ring(r, a, gg);
		if (r.empty()) {
			swap(g, gg);
			return true;
		}
		// next evaluation point
		xi = truncate1(xi*isqrt(isqrt(xi))*n73794, n27011);
	} while (--max_tries != 0);
	return false;
}

template<typename T> static bool
heur_gcd_z_priv(T& g, const T& a, const T& b, const unsigned max_tries = 66)
{
	typedef typename T::value_type ring_t;
	ring_t acont, bcont;
	T a_(a), b_(b);
	normalize_in_ring(a_, &acont);
	normalize_in_ring(b_, &bcont);
	const ring_t gc = gcd(acont, bcont);
	bool found = heur_gcd_z_pp(g, a_, b_, max_tries);
	if (found) {
		g *= gc;
		return true;
	}
	return false;
}

} // namespace GiNaC

#endif // ndef GINAC_UPOLY_HEUR_GCD_H
