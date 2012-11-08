/** @file exam_misc.cpp
 *
 *  Testing modular GCD.
 */

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

#include <cln/random.h>
#include <iostream>
#include <map>
#include <string>

#include "polynomial/upoly.h"
#include "polynomial/upoly_io.h"
#include "polynomial/mod_gcd.h"
#include "ginac.h"
using namespace GiNaC;

static upoly ex_to_upoly(const ex& e, const symbol& x);
static ex upoly_to_ex(const upoly& p, const symbol& x);

// make a univariate polynomial \in Z[x] of degree deg
static upoly make_random_upoly(const std::size_t deg);

static void run_test_once(const std::size_t deg)
{
	static const symbol xsym("x");

	const upoly a = make_random_upoly(deg);
	const upoly b = make_random_upoly(deg);

	upoly g;
	mod_gcd(g, a, b);

	ex ea = upoly_to_ex(a, xsym);
	ex eb = upoly_to_ex(b, xsym);
	ex eg = gcd(ea, eb);

	const upoly g_check = ex_to_upoly(eg, xsym);
	if (g != g_check) {
		std::cerr << "a = " << a << std::endl;
		std::cerr << "b = " << b << std::endl;
		std::cerr << "mod_gcd(a, b) = " << g << std::endl;
		std::cerr << "sr_gcd(a, b) = " << g_check << std::endl;
		throw std::logic_error("bug in mod_gcd");
	}
}

int main(int argc, char** argv)
{
	std::cout << "examining modular gcd. ";
	std::map<std::size_t, std::size_t> n_map;
	// run 256 tests with polynomials of degree 10
	n_map[10] = 256;
	// run 32 tests with polynomials of degree 100
	n_map[100] = 32;
	std::map<std::size_t, std::size_t>::const_iterator i = n_map.begin();
	for (; i != n_map.end(); ++i) {
		for (std::size_t k = 0; k < i->second; ++k)
			run_test_once(i->first);
	}
	return 0;
}

static upoly ex_to_upoly(const ex& e, const symbol& x)
{
	upoly p(e.degree(x) + 1);
	for (int i = 0; i <= e.degree(x); ++i)
		p[i] = cln::the<cln::cl_I>(ex_to<numeric>(e.coeff(x, i)).to_cl_N());
	return p;
}

static ex upoly_to_ex(const upoly& p, const symbol& x)
{
	exvector tv(p.size());
	for (std::size_t i = 0; i < p.size(); ++i)
		tv[i] = pow(x, i)*numeric(p[i]);
	return (new add(tv))->setflag(status_flags::dynallocated);
}

static upoly make_random_upoly(const std::size_t deg)
{
	static const cln::cl_I biggish("98765432109876543210");
	upoly p(deg + 1);
	for (std::size_t i = 0; i <= deg; ++i)
		p[i] = cln::random_I(biggish);

	// Make sure the leading coefficient is non-zero
	while (zerop(p[deg])) 
		p[deg] = cln::random_I(biggish);
	return p;
}
