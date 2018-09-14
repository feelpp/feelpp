/** @file remainder.h
 *
 *  Functions calculating remainders. */

/*
 *  GiNaC Copyright (C) 1999-2016 Johannes Gutenberg University Mainz, Germany
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

#ifndef GINAC_UPOLY_REMAINDER_H
#define GINAC_UPOLY_REMAINDER_H

#include "upoly.h"
#include "ring_traits.h"
#include "upoly_io.h"
#include "debug.h"

namespace GiNaC {

bool remainder_in_field(umodpoly& r, const umodpoly& a, const umodpoly& b);

/**
 * @brief Polynomial remainder for univariate polynomials over a ring. 
 *
 * Given two univariate polynomials \f$a, b \in R[x]\f$, where R is
 * a ring (presumably Z) computes the remainder @a r, which is
 * defined as \f$a = b q + r\f$. Returns true if the remainder is zero
 * and false otherwise.
 */
template<typename T>
bool remainder_in_ring(T& r, const T& a, const T& b)
{
	typedef typename T::value_type ring_t;
	if (degree(a) < degree(b)) {
		r = a;
		return false;
	}
	// N.B: don't bother to optimize division by constant

	r = a;
	const ring_t b_lcoeff = lcoeff(b);
	for (std::size_t k = a.size(); k-- >= b.size(); ) {

		// r -= r_k/b_n x^{k - n} b(x)
		if (zerop(r[k]))
			continue;

		const ring_t qk = truncate1(r[k], b_lcoeff);

		// Why C++ is so off-by-one prone?
		for (std::size_t j = k, i = b.size(); i-- != 0; --j) {
			if (zerop(b[i]))
				continue;
			r[j] = r[j] - qk*b[i];
		}

		if (!zerop(r[k])) {
			// division failed, don't bother to continue
			break;
		}
	}

	// Canonicalize the remainder: remove leading zeros. We can't say
	// anything about the degree of the remainder here.
	canonicalize(r);
	return r.empty();
}

} // namespace GiNaC

#endif // GINAC_UPOLY_REMAINDER_H
