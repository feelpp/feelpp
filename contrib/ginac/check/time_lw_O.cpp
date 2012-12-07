/** @file time_lw_O.cpp
 *
 *  Test O1 from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
 *  Lewis and Michael Wester. */

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
#include <vector>
using namespace std;

static const bool do_test2 = false;  // set to true in order to run this beast

static const symbol a1("a1"), a2("a2"), a3("a3"), a4("a4"), a5("a5"), a6("a6");
static const symbol b1("b1"), b2("b2"), b3("b3"), b4("b4"), b5("b5"), b6("b6");
static const symbol c1("c1"), c2("c2"), c3("c3"), c4("c4"), c5("c5"), c6("c6");

static const ex det1()
{
	matrix d1(15,15);
	d1 = a6, a5, a4, a3, a2, a1, 0,  0,  0,  0,  0,  0,  0,  0,  0,
	     0,  0,  a6, 0,  a5, a4, 0,  a3, a2, a1, 0,  0,  0,  0,  0,
	     0,  a6, 0,  a5, a4, 0,  a3, a2, a1, 0,  0,  0,  0,  0,  0,
	     0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1, 0,  0,
	     0,  0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1, 0,
	     0,  0,  0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1,
	     0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1, 0,  0,
	     0,  0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1, 0,
	     0,  b6, 0,  b5, b4, 0,  b3, b2, b1, 0,  0,  0,  0,  0,  0,
	     0,  0,  b6, 0,  b5, b4, 0,  b3, b2, b1, 0,  0,  0,  0,  0,
	     0,  0,  0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1,
	     0,  0,  0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1,
	     0,  0,  c6, 0,  c5, c4, 0,  c3, c2, c1, 0,  0,  0,  0,  0,
	     0,  c6, 0,  c5, c4, 0,  c3, c2, c1, 0,  0,  0,  0,  0,  0,
	     0,  0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1, 0;

	return d1.determinant();
}

static const ex det2()
{
	matrix d2(15,15);
	d2 = b6, b5, b4, b3, b2, b1, 0,  0,  0,  0,  0,  0,  0,  0,  0,
	     0,  0,  b6, 0,  b5, b4, 0,  b3, b2, b1, 0,  0,  0,  0,  0,
	     0,  b6, 0,  b5, b4, 0,  b3, b2, b1, 0,  0,  0,  0,  0,  0,
	     0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1, 0,  0,
	     0,  0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1, 0,
	     0,  0,  0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1,
	     0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1, 0,  0,
	     0,  0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1, 0,
	     0,  c6, 0,  c5, c4, 0,  c3, c2, c1, 0,  0,  0,  0,  0,  0,
	     0,  0,  c6, 0,  c5, c4, 0,  c3, c2, c1, 0,  0,  0,  0,  0,
	     0,  0,  0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1,
	     0,  0,  0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1,
	     0,  0,  a6, 0,  a5, a4, 0,  a3, a2, a1, 0,  0,  0,  0,  0,
	     0,  a6, 0,  a5, a4, 0,  a3, a2, a1, 0,  0,  0,  0,  0,  0,
	     0,  0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1, 0;

	return d2.determinant();
}

static const ex det3()
{
	matrix d3(15,15);
	d3 = c6, c5, c4, c3, c2, c1, 0,  0,  0,  0,  0,  0,  0,  0,  0,
	     0,  0,  c6, 0,  c5, c4, 0,  c3, c2, c1, 0,  0,  0,  0,  0,
	     0,  c6, 0,  c5, c4, 0,  c3, c2, c1, 0,  0,  0,  0,  0,  0,
	     0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1, 0,  0,
	     0,  0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1, 0,
	     0,  0,  0,  0,  0,  c6, 0,  0,  c5, c4, 0,  0,  c3, c2, c1,
	     0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1, 0,  0,
	     0,  0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1, 0,
	     0,  a6, 0,  a5, a4, 0,  a3, a2, a1, 0,  0,  0,  0,  0,  0,
	     0,  0,  a6, 0,  a5, a4, 0,  a3, a2, a1, 0,  0,  0,  0,  0,
	     0,  0,  0,  0,  0,  a6, 0,  0,  a5, a4, 0,  0,  a3, a2, a1,
	     0,  0,  0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1,
	     0,  0,  b6, 0,  b5, b4, 0,  b3, b2, b1, 0,  0,  0,  0,  0,
	     0,  b6, 0,  b5, b4, 0,  b3, b2, b1, 0,  0,  0,  0,  0,  0,
	     0,  0,  0,  0,  b6, 0,  0,  b5, b4, 0,  0,  b3, b2, b1, 0;

	return d3.determinant();
}

// The results of test_O1 will be needed for test_O2:
static ex d1, d2, d3;

static unsigned test_O1()
{
	d1 = det1();  cout << '.' << flush;
	d2 = det2();  cout << '.' << flush;
	d3 = det3();  cout << '.' << flush;
	unsigned nops1 = nops(d1);
	unsigned nops2 = nops(d2);
	unsigned nops3 = nops(d3);

	if ((nops1 != 37490) || (nops2 != 37490) || (nops3 != 37490)) {
		clog << "Determinants were miscalculated" << endl;
		return 1;
	}
	return 0;
}

static unsigned test_O2()
{
	const ex gcd1 = gcd( d1, d2 );  cout << '.' << flush;
	const ex resultant = gcd( gcd1, d3 );  cout << '.' << flush;
	if (nops(resultant) != 21894) {
		clog << "Resultant was miscalculated" << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_O()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;

	cout << "timing Lewis-Wester test O1 (three 15x15 dets)" << flush;

	rolex.start();
	// correct for very small times:
	do {
		result = test_O1();
		++count;
	} while ((time=rolex.read())<0.1 && !result);

	if (result)
		return result;
	
	cout << time/(3*count) << "s (average)" << endl;

	cout << "timing Lewis-Wester test O2 (Resultant)" << flush;

	if (do_test2) {
		rolex.reset();
		result += test_O2();

		cout << rolex.read() << "s (combined)" << endl;
	} else {
		cout << " disabled" << endl;
	}

	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_O();
}
