/** @file smod_helpers.h
 *
 *  Utility functions for modular arithmetic. */

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

#ifndef GINAC_POLYNOMIAL_SMOD_HELPERS_H
#define GINAC_POLYNOMIAL_SMOD_HELPERS_H

#include "ex.h"
#include "numeric.h"
#include "debug.h"

#include <cln/integer.h>
#include <cln/integer_io.h>

namespace GiNaC {

/// Z -> Z_p (in the symmetric representation)
static inline cln::cl_I smod(const cln::cl_I& a, long p)
{
	const cln::cl_I p2 = cln::cl_I(p >> 1);
	const cln::cl_I m = cln::mod(a, p);
	const cln::cl_I m_p = m - cln::cl_I(p);
	const cln::cl_I ret = m > p2 ? m_p : m;
	return ret;
}

static inline cln::cl_I recip(const cln::cl_I& a, long p_)
{
	cln::cl_I p(p_);
	cln::cl_I u, v;
	const cln::cl_I g = xgcd(a, p, &u, &v);
	cln::cl_I ret = smod(u, p_);
	cln::cl_I chck = smod(a*ret, p_);
	bug_on(chck != 1, "miscomputed recip(" << a << " (mod " << p_ << "))");
	return ret;

}

static inline numeric recip(const numeric& a_, long p)
{
	const cln::cl_I a = cln::the<cln::cl_I>(a_.to_cl_N());
	const cln::cl_I ret = recip(a, p);
	return numeric(ret);
}

static inline cln::cl_I to_cl_I(const ex& e)
{
	bug_on(!is_a<numeric>(e), "argument should be an integer");
	bug_on(!e.info(info_flags::integer),
		"argument should be an integer");
	return cln::the<cln::cl_I>(ex_to<numeric>(e).to_cl_N());
}

struct random_modint
{
	typedef long value_type;
	const value_type p;
	const value_type p_2;

	random_modint(const value_type& p_) : p(p_), p_2((p >> 1))
	{ }
	value_type operator()() const
	{
		do {
			cln::cl_I tmp_ = cln::random_I(p);
			value_type tmp = cln::cl_I_to_long(tmp_);
			if (tmp > p_2)
				tmp -= p;
			return tmp;
		} while (true);
	}

};

} // namespace GiNaC

#endif // GINAC_POLYNOMIAL_SMOD_HELPERS_H
