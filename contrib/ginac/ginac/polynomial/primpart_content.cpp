/** @file primpart_content.cpp
 *
 *  Function to find primitive part of a multivariate polynomial. */

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

#include "ex.h"
#include "numeric.h"
#include "collect_vargs.h"
#include "euclid_gcd_wrap.h"
#include "divide_in_z_p.h"
#include "debug.h"

namespace GiNaC {

/**
 * Compute the primitive part and the content of a modular multivariate
 * polynomial e \in Z_p[x_n][x_0, \ldots, x_{n-1}], i.e. e is considered
 * a polynomial in variables x_0, \ldots, x_{n-1} with coefficients being
 * modular polynomials Z_p[x_n]
 * @param e polynomial to operate on
 * @param pp primitive part of @a e, will be computed by this function
 * @param c content (in the sense described above) of @a e, will be
 *        computed by this function
 * @param vars variables x_0, \ldots, x_{n-1}, x_n
 * @param p modulus
 */
void primpart_content(ex& pp, ex& c, ex e, const exvector& vars,
		      const long p)
{
	static const ex ex1(1);
	static const ex ex0(0);
	e = e.expand();
	if (e.is_zero()) {
		pp = ex0;
		c = ex1;
		return;
	}
	exvector rest_vars = vars;
	rest_vars.pop_back();
	ex_collect_t ec;
	collect_vargs(ec, e, rest_vars);

	if (ec.size() == 1) {
		// the input polynomial factorizes into 
		// p_1(x_n) p_2(x_0, \ldots, x_{n-1})
		c = ec.rbegin()->second;
		ec.rbegin()->second = ex1;
		pp = ex_collect_to_ex(ec, rest_vars).expand().smod(numeric(p));
		return;
	}

	// Start from the leading coefficient (which is stored as a last
	// element of the terms array)
	ex_collect_t::reverse_iterator i = ec.rbegin();
	ex g = i->second;
	// there are at least two terms, so it's safe to...
	++i;
	while (i != ec.rend() && !g.is_equal(ex1)) {
		g = euclid_gcd(i->second, g, vars.back(), p);
		++i;
	}
	if (g.is_equal(ex1)) {
		pp = e;
		c = ex1;
		return;
	}
	exvector mainvar;
	mainvar.push_back(vars.back());
	for (i = ec.rbegin(); i != ec.rend(); ++i) {
		ex tmp(0);
		if (!divide_in_z_p(i->second, g, tmp, mainvar, p))
			throw std::logic_error(std::string(__func__) +
					": bogus division failure");
		i->second = tmp;
	}

	pp = ex_collect_to_ex(ec, rest_vars).expand().smod(numeric(p));
	c = g;
}

} // namespace GiNaC
