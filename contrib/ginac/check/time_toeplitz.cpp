/** @file time_toeplitz.cpp
 *
 *  Calculates determinants of dense symbolic Toeplitz materices.
 *  For 4x4 our matrix would look like this:
 *  [[a,b,a+b,a^2+a*b+b^2], [b,a,b,a+b], [a+b,b,a,b], [a^2+a*b+b^2,a+b,b,a]]
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
#include <vector>
using namespace std;

static unsigned toeplitz_det(unsigned size)
{
	unsigned result = 0;
	const symbol a("a"), b("b");
	ex p[10] = {ex("a",lst(a,b)),
	           ex("b",lst(a,b)),
	           ex("a+b",lst(a,b)),
	           ex("a^2+a*b+b^2",lst(a,b)),
	           ex("a^3+a^2*b-a*b^2+b^3",lst(a,b)),
	           ex("a^4+a^3*b+a^2*b^2+a*b^3+b^4",lst(a,b)),
	           ex("a^5+a^4*b+a^3*b^2-a^2*b^3+a*b^4+b^5",lst(a,b)),
	           ex("a^6+a^5*b+a^4*b^2+a^3*b^3+a^2*b^4+a*b^5+b^6",lst(a,b)),
	           ex("a^7+a^6*b+a^5*b^2+a^4*b^3-a^3*b^4+a^2*b^5+a*b^6+b^7",lst(a,b)),
	           ex("a^8+a^7*b+a^6*b^2+a^5*b^3+a^4*b^4+a^3*b^5+a^2*b^6+a*b^7+b^8",lst(a,b))
	};

	// construct Toeplitz matrix (diagonal structure: [[x,y,z],[y,x,y],[z,y,x]]):
	matrix M(size,size);
	for (unsigned ro=0; ro<size; ++ro) {
		for (unsigned nd=ro; nd<size; ++nd) {
			M.set(nd-ro,nd,p[ro]);
			M.set(nd,nd-ro,p[ro]);
		}
	}

	// compute determinant:
	ex tdet = M.determinant();

	// dirty consistency check of result:
	if (!tdet.subs(a==0).subs(b==0).is_zero()) {
		clog << "Determaint of Toeplitz matrix " << endl
		     << "M==" << M << endl
		     << "was miscalculated: det(M)==" << tdet << endl;
		++result;
	}

	return result;
}

unsigned time_toeplitz()
{
	unsigned result = 0;

	cout << "timing determinant of polyvariate symbolic Toeplitz matrices" << flush;

	vector<unsigned> sizes;
	vector<double> times;
	timer longines;

	sizes.push_back(7);
	sizes.push_back(8);
	sizes.push_back(9);
	sizes.push_back(10);

	for (vector<unsigned>::iterator i=sizes.begin(); i!=sizes.end(); ++i) {
		int count = 1;
		longines.start();
		result += toeplitz_det(*i);
		// correct for very small times:
		while (longines.read()<0.1) {
			toeplitz_det(*i);
			++count;
		}
		times.push_back(longines.read()/count);
		cout << '.' << flush;
	}

	// print the report:
	cout << endl << "	dim:   ";
	for (vector<unsigned>::iterator i=sizes.begin(); i!=sizes.end(); ++i)
		cout << '\t' << *i << 'x' << *i;
	cout << endl << "	time/s:";
	for (vector<double>::iterator i=times.begin(); i!=times.end(); ++i)
		cout << '\t' << *i;
	cout << endl;

	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_toeplitz();
}
