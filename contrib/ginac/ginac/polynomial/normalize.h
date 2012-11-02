/** @file normalize.h
 *
 *  Functions to normalize polynomials in a field. */

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

#ifndef GINAC_UPOLY_NORMALIZE_H
#define GINAC_UPOLY_NORMALIZE_H

#include "upoly.h"
#include "ring_traits.h"
#include "debug.h"

namespace GiNaC {

/// Make the univariate polynomial @a a \in F[x] unit normal.
/// F should be a field.
/// Returns true if the polynomial @x is already unit normal, and false
/// otherwise.
static bool normalize_in_field(umodpoly& a, cln::cl_MI* content_ = 0)
{
	if (a.size() == 0)
		return true;
	if (lcoeff(a) == the_one(a[0])) {
		if (content_)
			*content_ = the_one(a[0]);
		return true;
	}

	const cln::cl_MI lc_1 = recip(lcoeff(a));
	for (std::size_t k = a.size(); k-- != 0; )
		a[k] = a[k]*lc_1;
	if (content_)
		*content_ = lc_1;
	return false;
}

/// Make the univariate polynomial @a x unit normal. This version is used
/// for rings which are not fields. 
/// Returns true if the polynomial @x is already unit normal, and false
/// otherwise.
template<typename T> bool
normalize_in_ring(T& x, typename T::value_type* content_ = 0, int* unit_ = 0)
{
	typedef typename T::value_type ring_t;
	static const ring_t one(1);
	if (x.empty())
		return true;

	bool something_changed = false;
	if (minusp(lcoeff(x))) {
		something_changed = true;
		if (unit_)
			*unit_ = -1;
		for (std::size_t i = x.size(); i-- != 0; )
			x[i] = - x[i];
	}

	if (degree(x) == 0) {
		if (content_)
			*content_ = x[0];
		if (x[0] == one)
			return something_changed;
		x[0] = one;
		return false; // initial polynomial was unit normal
	}
		
	// Compute the gcd of coefficients
	ring_t content = lcoeff(x);
	// We want this function to be fast when applied to unit normal
	// polynomials. Hence we start from the leading coefficient.
	for (std::size_t i = x.size() - 1; i-- != 0; ) {
		if (content == one) {
			if (content_)
				*content_ = one;
			return something_changed;
		}
		content = gcd(x[i], content);
	}

	if (content == one) {
		if (content_)
			*content_ = one;
		return something_changed;
	}
	
	for (std::size_t i = x.size(); i-- != 0; )
		x[i] = exquo(x[i], content);

	if (content_)
		*content_ = content;
	return false; // initial polynomial was not unit normal
}

} // namespace GiNaC
	
#endif // GINAC_UPOLY_NORMALIZE_H
