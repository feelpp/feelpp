/** @file gcd_euclid.h
 *
 *  GCD using Euclidean algorithm. */

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

#ifndef GINAC_GCD_EUCLID_H
#define GINAC_GCD_EUCLID_H

#include "upoly.h"
#include "remainder.h"
#include "normalize.h"
#include "debug.h"
#include "upoly_io.h"

namespace GiNaC {

static bool
gcd_euclid(umodpoly& c, umodpoly /* passed by value */ a, umodpoly b)
{
	if (a.size() == 0) {
		c.clear();
		return true;
	}
	if (b.size() == 0) {
		c.clear();
		return true;
	}
	bug_on(a[0].ring()->modulus != b[0].ring()->modulus,
		"different moduli");

	normalize_in_field(a);
	normalize_in_field(b);
	if (degree(a) < degree(b))
		std::swap(a, b);

	umodpoly r;
	while (b.size() != 0) {
		remainder_in_field(r, a, b); 
		a = b;
		b = r;
	}
	normalize_in_field(a);
	c = a;
	return false;
}

} // namespace GiNaC

#endif // ndef GINAC_GCD_EUCLID_H
