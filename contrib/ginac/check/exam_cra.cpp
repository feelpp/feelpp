/** @file exam_cra.cpp
 *
 *  Test of Chinese remainder algorithm. */

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

#include "polynomial/cra_garner.h"

#include <cln/integer.h>
#include <cln/integer_io.h>
#include <cln/random.h>
#include <cln/numtheory.h>
using namespace cln;
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;

/// Generate a sequences of primes p_i such that \prod_i p_i < limit
static std::vector<cln::cl_I>
make_random_moduli(const cln::cl_I& limit);

static std::vector<cln::cl_I>
calc_residues(const cln::cl_I& x, const std::vector<cln::cl_I>& moduli);

static void dump(const std::vector<cln::cl_I>& v);

/// Make @a n random relatively prime moduli, each < limit, make a 
/// random number x < \prod_{i=0}{n-1}, calculate residues, and 
/// compute x' by chinese remainder algorithm. Check if the result
/// of computation matches the original value x.
static void run_test_once(const cln::cl_I& lim)
{
	std::vector<cln::cl_I> moduli = make_random_moduli(lim);
	cln::cl_I x = random_I(lim) + 1;

	if (x > (lim >> 1))
		x = x - lim;

	std::vector<cln::cl_I> residues = calc_residues(x, moduli);
	cln::cl_I x_test;

	bool error = false;
	try {
		x_test = integer_cra(residues, moduli);
	} catch (std::exception& oops) {
		std::cerr << "Oops: " << oops.what() << std::endl;
		error = true;
	}

	if (x != x_test)
		error = true;

	if (error) {
		std::cerr << "Expected x = " << x << ", got " <<
			x_test << " instead" << std::endl;
		std::cerr << "moduli = ";
		dump(moduli);
		std::cerr << std::endl;
		std::cerr << "residues = ";
		dump(residues);
		std::cerr << std::endl;
		throw std::logic_error("bug in integer_cra?");
	}
}

static void run_test(const cln::cl_I& limit, const std::size_t ntimes)
{
	for (std::size_t i = 0; i < ntimes; ++i)
		run_test_once(limit);
}

int main(int argc, char** argv)
{
	typedef std::map<cln::cl_I, std::size_t> map_t;
	map_t the_map;
	// Run 1024 tests with native 32-bit numbers
	the_map[cln::cl_I(std::numeric_limits<int>::max())] = 1024;

	// Run 512 tests with native 64-bit integers
	if (sizeof(long) > sizeof(int)) 
		the_map[cln::cl_I(std::numeric_limits<long>::max())] = 512;

	// Run 32 tests with a bit bigger numbers
	the_map[cln::cl_I("987654321098765432109876543210")] = 32;

	std::cout << "examining Garner's integer chinese remainder algorithm " << std::flush;

	for (map_t::const_iterator i = the_map.begin(); i != the_map.end(); ++i)
		run_test(i->first, i->second);

	return 0;
}

static std::vector<cln::cl_I>
calc_residues(const cln::cl_I& x, const std::vector<cln::cl_I>& moduli)
{
	std::vector<cln::cl_I> residues(moduli.size());
	for (std::size_t i = moduli.size(); i-- != 0; )
		residues[i] = mod(x, moduli[i]);
	return residues;
}

static std::vector<cln::cl_I>
make_random_moduli(const cln::cl_I& limit)
{
	std::vector<cln::cl_I> moduli;
	cln::cl_I prod(1);
	cln::cl_I next = random_I(std::min(limit >> 1, cln::cl_I(128)));
	unsigned count = 0;
	do {
		cln::cl_I tmp = nextprobprime(next);
		next = tmp + random_I(cln::cl_I(10)) + 1;
		prod = prod*tmp;
		moduli.push_back(tmp);
		++count;
	} while (prod < limit || (count < 2));
	return moduli;
}

static void dump(const std::vector<cln::cl_I>& v)
{
	std::cerr << "[ ";
	for (std::size_t i = 0; i < v.size(); ++i)
		std::cerr << v[i] << " ";
	std::cerr << "]";
}
