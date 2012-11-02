/** @file time_parser.cpp
 *
 *  Time the parser. */

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

#include "ginac.h"
#include "timer.h"
using namespace GiNaC;

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

extern void randomify_symbol_serials();

/// make a string "1+x+2*x^2+...+n*x^n"
static string prepare_str(const unsigned n, const char x = 'x')
{
	ostringstream s;
	s << x;
	for (unsigned i = 2; i < n; i++)
		s << '+' << i << '*' << x << '^' << i; 
	return s.str();
}

static double benchmark_and_cmp(const string& srep)
{
	parser the_parser;
	timer RSD10;
	RSD10.start();
	ex e = the_parser(srep);
	const double t = RSD10.read();
	return t;
}

int main(int argc, char** argv)
{
	cout << "timing GiNaC parser..." << flush;
	randomify_symbol_serials();
	unsigned n_min = 1024;
	unsigned n_max = 32768;
	if (argc > 1)
		n_max = atoi(argv[1]);

	vector<double> times;
	vector<unsigned> ns;
	for (unsigned n = n_min; n <= n_max; n = n << 1) {
		string srep = prepare_str(n);
		const double t = benchmark_and_cmp(srep);
		times.push_back(t);
		ns.push_back(n);
	}

	cout << "OK" << endl;
	cout << "# terms  time, s" << endl;
	for (size_t i = 0; i < times.size(); i++)
		cout << " " << ns[i] << '\t' << times[i] << endl;
	return 0;
}
