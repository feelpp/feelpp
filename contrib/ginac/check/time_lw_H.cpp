/** @file time_lw_H.cpp
 *
 *  Test H from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
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

static unsigned test(unsigned n)
{
	matrix hilbert(n,n);

	for (unsigned r=0; r<n; ++r)
		for (unsigned c=0; c<n; ++c)
			hilbert.set(r,c,numeric(1,r+c+1));
	ex det = hilbert.determinant();

	/*
	   The closed form of the determinant of n x n Hilbert matrices is:
	
             n-1   /                   \
            ----- | pow(factorial(r),3) |
             | |  | ------------------- |
             | |  |    factorial(r+n)   |
            r = 0  \                   /
	*/

	ex hilbdet = 1;
	for (unsigned r=0; r<n; ++r)
		hilbdet *= pow(factorial(r),3)/(factorial(r+n));

	if (det != hilbdet) {
		clog << "determinant of " << n << "x" << n << " erroneously returned " << det << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_H()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;

	cout << "timing Lewis-Wester test H (det of 80x80 Hilbert)" << flush;

	rolex.start();
	// correct for very small times:
	do {
		result = test(80);
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
	return time_lw_H();
}
