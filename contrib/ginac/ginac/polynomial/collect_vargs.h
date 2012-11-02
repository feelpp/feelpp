/** @file collect_vargs.h
 *
 *  Interface to utility functions. */

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

#ifndef GINAC_COLLECT_VARGS_H
#define GINAC_COLLECT_VARGS_H

#include "ex.h"

#include <cln/integer.h>
#include <utility> // for std::pair
#include <vector>
#include <algorithm> // std::lexicographical_compare

namespace GiNaC {

typedef std::vector<int> exp_vector_t;

static inline bool operator<(const exp_vector_t& v1, const exp_vector_t& v2)
{
	return std::lexicographical_compare(v1.rbegin(), v1.rend(),
					    v2.rbegin(), v2.rend());
}

static inline bool operator>(const exp_vector_t& v1, const exp_vector_t& v2)
{
	if (v1 == v2)
		return false;
	return !(v1 < v2);
}

static inline bool zerop(const exp_vector_t& v)
{
	for (exp_vector_t::const_reverse_iterator i = v.rbegin(); i != v.rend(); ++i) {
		if (*i != 0) 
			return false;
	}
	return true;
}

typedef std::vector<std::pair<exp_vector_t, ex> > ex_collect_t;

extern void
collect_vargs(ex_collect_t& ec, const ex& e, const exvector& x);
extern ex
ex_collect_to_ex(const ex_collect_t& ec, const exvector& x);

/**
 * Leading coefficient of a multivariate polynomial e, considering it
 * as a multivariate polynomial in x_0, \ldots x_{n-1} with coefficients
 * being univariate polynomials in R[x_n] (where R is some ring)
 */
extern ex lcoeff_wrt(ex e, const exvector& x);


/**
 * Degree vector of a leading term of a multivariate polynomial.
 * (generalization of degree(expr, var))
 */
extern exp_vector_t degree_vector(ex e, const exvector& vars);

/**
 * Leading coefficient c \in R (where R = Z or Z_p) of a multivariate
 * polynomial e \in R[x_0, \ldots, x_n]
 */
extern cln::cl_I integer_lcoeff(const ex& e, const exvector& vars);

} // namespace GiNaC

#endif // ndef GINAC_COLLECT_VARGS_H
