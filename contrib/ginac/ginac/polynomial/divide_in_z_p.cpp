/** @file divide_in_z_p.cpp
 *
 *  Implementation of polynomial division in Z/Zp. */

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

#include "add.h"
#include "operators.h"
#include "power.h"
#include "smod_helpers.h"

namespace GiNaC {

/** 
 * Exact polynomial division of a, b \in Z_p[x_0, \ldots, x_n]
 * It doesn't check whether the inputs are proper polynomials, so be careful
 * of what you pass in.
 *  
 * @param a  first multivariate polynomial (dividend)
 * @param b  second multivariate polynomial (divisor)
 * @param q  quotient (returned)
 * @param var variables X iterator to first element of vector of symbols
 *
 * @return "true" when exact division succeeds (the quotient is returned in
 *          q), "false" otherwise.
 * @note @a p = 0 means the base ring is Z
 */
bool divide_in_z_p(const ex &a, const ex &b, ex &q, const exvector& vars, const long p)
{
	static const ex _ex1(1);
	if (b.is_zero())
		throw(std::overflow_error("divide_in_z: division by zero"));
	if (b.is_equal(_ex1)) {
		q = a;
		return true;
	}
	if (is_exactly_a<numeric>(a)) {
		if (is_exactly_a<numeric>(b)) {
			// p == 0 means division in Z
			if (p == 0) {
				const numeric tmp = ex_to<numeric>(a/b);
				if (tmp.is_integer()) {
					q = tmp;
					return true;
				} else
					return false;
			} else {
				q = (a*recip(ex_to<numeric>(b), p)).smod(numeric(p));
				return true;
			}
		} else
			return false;
	}
	if (a.is_equal(b)) {
		q = _ex1;
		return true;
	}

	// Main symbol
	const ex &x = vars.back();

	// Compare degrees
	int adeg = a.degree(x), bdeg = b.degree(x);
	if (bdeg > adeg)
		return false;

	// Polynomial long division (recursive)
	ex r = a.expand();
	if (r.is_zero())
		return true;
	int rdeg = adeg;
	ex eb = b.expand();
	ex blcoeff = eb.coeff(x, bdeg);
	exvector v;
	v.reserve(std::max(rdeg - bdeg + 1, 0));
	exvector rest_vars(vars);
	rest_vars.pop_back();
	while (rdeg >= bdeg) {
		ex term, rcoeff = r.coeff(x, rdeg);
		if (!divide_in_z_p(rcoeff, blcoeff, term, rest_vars, p))
			break;
		term = (term*power(x, rdeg - bdeg)).expand();
		v.push_back(term);
		r = (r - term*eb).expand();
		if (p != 0)
			r = r.smod(numeric(p));
		if (r.is_zero()) {
			q = (new add(v))->setflag(status_flags::dynallocated);
			return true;
		}
		rdeg = r.degree(x);
	}
	return false;
}

} // namespace GiNaC
