/** @file collect_vargs.cpp
 *
 *  Utility functions. */

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
#include "mul.h"
#include "operators.h"
#include "power.h"
#include "collect_vargs.h"
#include "smod_helpers.h"
#include "debug.h"

#include <algorithm>
#include <cln/integer.h>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>

namespace GiNaC {

typedef std::map<exp_vector_t, ex> ex_collect_priv_t;

static void 
collect_vargs(ex_collect_priv_t& ec, ex e, const exvector& vars);
static void
collect_term(ex_collect_priv_t& ec, const ex& e, const exvector& vars);
static void wipe_out_zeros(ex_collect_priv_t& ec);

template<typename T, typename CoeffCMP>
struct compare_terms
{
	const CoeffCMP& coeff_cmp;
	explicit compare_terms(const CoeffCMP& coeff_cmp_) : coeff_cmp(coeff_cmp_)
	{ }
	inline bool operator()(const T& t1, const T& t2) const
	{
		bool exponent_is_less =
			std::lexicographical_compare(t1.first.rbegin(),
						     t1.first.rend(),
						     t2.first.rbegin(),
						     t2.first.rend());
		if (exponent_is_less)
			return true;

		if ((t1.first == t2.first) &&
				coeff_cmp(t2.second, t2.second))
			return true;
		return false;
	}
};

template<typename T, typename CoeffCMP>
static compare_terms<T, CoeffCMP>
make_compare_terms(const T& dummy, const CoeffCMP& coeff_cmp)
{
	return compare_terms<T, CoeffCMP>(coeff_cmp);
}

void collect_vargs(ex_collect_t& ec, const ex& e, const exvector& vars)
{
	ex_collect_priv_t ecp;
	collect_vargs(ecp, e, vars);
	ec.reserve(ecp.size());
	std::copy(ecp.begin(), ecp.end(), std::back_inserter(ec));
	std::sort(ec.begin(), ec.end(),
		  make_compare_terms(*ec.begin(), ex_is_less()));
}

static void 
collect_vargs(ex_collect_priv_t& ec, ex e, const exvector& vars)
{
	e = e.expand();
	if (e.is_zero()) {
		ec.clear();
		return;
	}

	if (!is_a<add>(e)) {
		collect_term(ec, e, vars);
		return;
	}

	for (const_iterator i = e.begin(); i != e.end(); ++i)
		collect_term(ec, *i, vars);

	wipe_out_zeros(ec);
}

static void
collect_term(ex_collect_priv_t& ec, const ex& e, const exvector& vars)
{
	if (e.is_zero())
		return;
	static const ex ex1(1);
	exp_vector_t key(vars.size());
	ex pre_coeff = e;
	for (std::size_t i = 0; i < vars.size(); ++i) {
		const int var_i_pow = pre_coeff.degree(vars[i]);
		key[i] = var_i_pow;
		pre_coeff = pre_coeff.coeff(vars[i], var_i_pow);
	}
	ex_collect_priv_t::iterator i = ec.find(key);
	if (i != ec.end())
		i->second += pre_coeff;
	else
		ec.insert(ex_collect_priv_t::value_type(key, pre_coeff));
}

static void wipe_out_zeros(ex_collect_priv_t& m)
{
	ex_collect_priv_t::iterator i = m.begin();
	while (i != m.end()) {
		// be careful to not invalide iterator, use post-increment
		// for that, see e.g.
		// http://coding.derkeiler.com/Archive/C_CPP/comp.lang.cpp/2004-02/0502.html
		if (i->second.is_zero())
			m.erase(i++);
		else
			++i;
	}
}

GiNaC::ex
ex_collect_to_ex(const ex_collect_t& ec, const GiNaC::exvector& vars)
{
	exvector ev;
	ev.reserve(ec.size());
	for (std::size_t i = 0; i < ec.size(); ++i) {
		exvector tv;
		tv.reserve(vars.size() + 1);
		for (std::size_t j = 0; j < vars.size(); ++j) {
			const exp_vector_t& exp_vector(ec[i].first);

			bug_on(exp_vector.size() != vars.size(),
				"expected " << vars.size() << " variables, "
				"expression has " << exp_vector.size() << " instead");

			if (exp_vector[j] != 0)
				tv.push_back(power(vars[j], exp_vector[j]));
		}
		tv.push_back(ec[i].second);
		ex tmp = (new mul(tv))->setflag(status_flags::dynallocated);
		ev.push_back(tmp);
	}
	ex ret = (new add(ev))->setflag(status_flags::dynallocated);
	return ret;
}

ex lcoeff_wrt(ex e, const exvector& x)
{
	static const ex ex0(0);
	e = e.expand();
	if (e.is_zero())
		return ex0;

	ex_collect_t ec;
	collect_vargs(ec, e, x);
	return ec.rbegin()->second;
}

exp_vector_t degree_vector(ex e, const exvector& vars)
{
	e = e.expand();
	exp_vector_t dvec(vars.size());
	for (std::size_t i = vars.size(); i-- != 0; ) {
		const int deg_i = e.degree(vars[i]);
		e = e.coeff(vars[i], deg_i);
		dvec[i] = deg_i;
	}
	return dvec;
}

cln::cl_I integer_lcoeff(const ex& e, const exvector& vars)
{
	ex_collect_t ec;
	collect_vargs(ec, e, vars);
	if (ec.size() == 0)
		return cln::cl_I(0);
	ex lc = ec.rbegin()->second;
	bug_on(!is_a<numeric>(lc), "leading coefficient is not an integer");
	bug_on(!lc.info(info_flags::integer),
		"leading coefficient is not an integer");

	return to_cl_I(lc);
}

} // namespace GiNaC
