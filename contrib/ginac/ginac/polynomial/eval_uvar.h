/** @file eval_uvar.h
 *
 *  Numerical evaluation of univariate polynomials. */

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

#ifndef GINAC_EVAL_UPOLY_H
#define GINAC_EVAL_UPOLY_H

#include "upoly.h"
#include "ring_traits.h"

namespace GiNaC {

/// Evaluate the polynomial using Horner rule.
/// TODO: 
/// - a better algorithm for small polynomials (use SIMD instructions)
/// - a better algorithm for large polynomials (use Karatsuba trick)
/// - a better algorithm for modular polynomials (especially for small
///     moduli and GFN)
template<typename T> static typename T::value_type
eval(const T& p, const typename T::value_type& x)
{
	// p(x) = c_n x^n + c_{n-1} x^{n-1} + \ldots + c_0 =
	// c_0 + x (c_1 + x (c_2 + x ( \ldots (c_{n-1} + c_n x) \ldots )))
	// (AKA Horner rule)
	typedef typename T::value_type ring_t;
	if (p.empty())
		return 0;
	if (degree(p) == 0)
		return p[0];

	ring_t y = lcoeff(p);
	// read the formula above from the right to the left
	for (std::size_t i = p.size() - 1; i-- != 0; )
		y = x*y + p[i];

	return y;
}

} // namespace GiNaC

#endif // ndef GINAC_EVAL_UPOLY_H
