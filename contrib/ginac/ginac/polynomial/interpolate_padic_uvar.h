/** @file interpolate_padic_uvar.h
 *
 *  Utility function for interpolation. */

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

#ifndef GINAC_UPOLY_INTERPOLATE_PADIC_H
#define GINAC_UPOLY_INTERPOLATE_PADIC_H

#include "ring_traits.h"

namespace cln {

static inline cln::cl_I smod(const cln::cl_I& x, const cln::cl_I y)
{
	cln::cl_I r = mod(x, y);
	const cln::cl_I y_2 = y >> 1;
	if (r > y_2)
		r = r - y;
	return r;
}

} // namespace cln

namespace GiNaC {

template<typename T> static void
interpolate(T& g, const typename T::value_type& gamma,
	    const typename T::value_type& modulus,
	    const std::size_t degree_hint = 1)
{
	typedef typename T::value_type ring_t;
	g.clear();
	g.reserve(degree_hint);
	ring_t e = gamma;
	while (!zerop(e)) {
		const ring_t gi = smod(e, modulus);
		g.push_back(gi);
		e = exquo(e - gi, modulus);
	}
}

} // namespace GiNaC

#endif // ndef GINAC_UPOLY_INTERPOLATE_PADIC_H
