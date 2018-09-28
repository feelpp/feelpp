/** @file normalize.h
 *
 *  Functions to normalize polynomials in a field. */

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

#include "normalize.h"

namespace GiNaC {

/// Make the univariate polynomial @a a \in F[x] unit normal.
/// F should be a field.
/// Returns true if the polynomial @x is already unit normal, and false
/// otherwise.
bool normalize_in_field(umodpoly& a, cln::cl_MI* content_)
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

} // namespace GiNaC
