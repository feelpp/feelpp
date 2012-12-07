/** @file ring_traits.h
 *
 *  Functions for polynomial ring arithmetic. */

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

#ifndef GINAC_RING_TRAITS_H
#define GINAC_RING_TRAITS_H

#include <cln/integer.h>
#include <cln/modinteger.h>

namespace cln {

static inline cln::cl_I div(const cln::cl_I& x, const cln::cl_I& y)
{
	return cln::exquo(x, y);
}

/// Exact integer division.
/// Check if y divides x, if yes put the quotient into q, otherwise don't
/// touch q. Returns true if y divides x and false if not.
static inline bool div(cln::cl_I& q, const cln::cl_I& x, const cln::cl_I& y)
{
	const cln::cl_I_div_t qr = cln::truncate2(x, y);
	if (zerop(qr.remainder)) {
		q = qr.quotient;
		return true;
	}
	return false;
}

static inline cln::cl_I get_ring_elt(const cln::cl_I& sample, const int val)
{
	return cln::cl_I(val);
}

static inline cln::cl_MI get_ring_elt(const cln::cl_MI& sample, const int val)
{
	return sample.ring()->canonhom(val);
}

template<typename T>
static inline T the_one(const T& sample)
{
	return get_ring_elt(sample, 1);
}

} // namespace cln

#endif // GINAC_RING_TRAITS_H
