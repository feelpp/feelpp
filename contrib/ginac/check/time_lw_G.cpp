/** @file time_lw_G.cpp
 *
 *  Test G from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
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

static unsigned test()
{
	symbol x("x");
	symbol y("y");
	symbol z("z");
	
	ex p = expand(pow(7*y*pow(x*z,2)-3*x*y*z+11*(x+1)*pow(y,2)+5*z+1,4)
	              *pow(3*x-7*y+2*z-3,5));
	ex q = expand(pow(7*y*pow(x*z,2)-3*x*y*z+11*(x+1)*pow(y,2)+5*z+1,3)
	              *pow(3*x-7*y+2*z+3,6));
	ex result = gcd(p,q);
	if (result.expand()!=expand(pow(7*y*pow(x*z,2)-3*x*y*z+11*(x+1)*pow(y,2)+5*z+1,3))) {
		clog << "gcd(expand((7*y*x^2*z^2-3*x*y*z+11*(x+1)*y^2+5*z+1)^4*(3*x-7*y+2*z-3)^5),expand((7*y*x^2*z^2-3*x*y*z+11*(x+1)*y^2+5*z+1)^3*(3*x-7*y+2*z+3)^6)) erroneously returned " << result << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_G()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;
	
	cout << "timing Lewis-Wester test G (gcd of 3-var polys)" << flush;
	
	rolex.start();
	// correct for very small times:
	do {
		result = test();
		++count;
	} while ((time=rolex.read())<0.1 && !result);
	cout << '.' << flush;
	cout << time/count << 's' << endl;
	
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_G();
}
