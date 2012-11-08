/** @file heur_gcd_bug.cpp
 *
 *  heur_gcd_oops.cpp Check for a bug in heur_gcd().
 *
 *  heur_gcd() did not check if the arguments are integer polynomials
 *  (and did not convert them to integer polynomials), which lead to
 *  endless loop or (even worse) wrong result. */

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
using namespace GiNaC;

#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
	cout << "checking if heur_gcd() can cope with rational polynomials. ";
	const symbol x("x");
	const ex _ex1(1);
	ex a1 = x + numeric(5, 4);
	ex a2 = x + numeric(5, 2);
	ex b =  pow(x, 2) + numeric(15, 4)*x + numeric(25, 8);
	// note: both a1 and a2 divide b
	
	// a2 divides b, so cofactor of a2 should be a (rational) number
	ex ca2, cb2;
	ex g2 = gcd(a2, b, &ca2, &cb2);
	if (!is_a<numeric>(ca2)) {
		cerr << "gcd(" << a2 << ", " << b << ") was miscomputed" << endl;
		return 1;
	}
	ex ca1, cb1;
	// a1 divides b, so cofactor of a1 should be a (rational) number
	ex g1 = gcd(a1, b, &ca1, &cb1);
	if (!is_a<numeric>(ca1)) {
		cerr << "gcd(" << a1 << ", " << b << ") was miscomputed" << endl;
		return 1;
	}
	return 0;
}
