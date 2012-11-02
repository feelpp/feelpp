/** @file time_fateman_expand.cpp
 *
 *  Time for polynomial expansion of (x+y+z+1)^20 * ((x+y+z+1)^20+1).
 *  This test was suggested by Richard J. Fateman as a benchmark for programs
 *  to multiply sparse polynomials fast.
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

#include "ginac.h"
#include "timer.h"
using namespace GiNaC;

#include <iostream>
using namespace std;

static unsigned test()
{
	unsigned result = 0;
	const symbol x("x"), y("y"), z("z");

	const ex p = pow(x+y+z+1, 20);

	const ex hugesum = expand(p * (p+1));

	if (hugesum.nops()!=12341) {
		clog << "(x+y+z+1)^20 * ((x+y+z+1)^20+1) was miscomputed!" << endl;
		++result;
	}

	return result;
}

unsigned time_fateman_expand()
{
	unsigned result = 0;
	unsigned count = 0;
	timer concord;
	double time = .0;

	cout << "timing Fateman's polynomial expand benchmark" << flush;

	concord.start();
	// correct for very small times:
	do {
		result = test();
		++count;
	} while ((time=concord.read())<0.1 && !result);
	cout << '.' << flush;

	cout << time/count << 's' << endl;

	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_fateman_expand();
}
