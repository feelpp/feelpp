/** @file time_vandermonde.cpp
 *
 *  Calculates determinants of dense symbolic Vandermonde materices with
 *  monomials in one single variable as entries.
 *  For 4x4 our matrix would look like this:
 *  [[1,a,a^2,a^3], [1,-a,a^2,-a^3], [1,a^2,a^4,a^6], [1,-a^2,a^4,-a^6]]
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

static unsigned vandermonde_det(unsigned size)
{
	unsigned result = 0;
	const symbol a("a");

	// construct Vandermonde matrix:
	matrix M(size,size);
	for (unsigned ro=0; ro<size; ++ro) {
		for (unsigned co=0; co<size; ++co) {
			if (ro%2)
				M(ro,co) = pow(-pow(a,1+ro/2),co);
			else
				M(ro,co) = pow(pow(a,1+ro/2),co);
		}
	}

	// compute determinant:
	ex det = M.determinant();

	// check the result:
	ex vanddet = 1;
	for (unsigned i=0; i<size; ++i)
		for (unsigned j=0; j<i; ++j)
			vanddet *= M(i,1) - M(j,1);

	if (expand(det - vanddet) != 0) {
		clog << "Determaint of Vandermonde matrix " << endl
		     << "M==" << M << endl
		     << "was miscalculated: det(M)==" << det << endl;
		++result;
	}

	return result;
}

unsigned time_vandermonde()
{
	unsigned result = 0;
	
	cout << "timing determinant of univariate symbolic Vandermonde matrices" << flush;
	
	vector<unsigned> sizes;
	vector<double> times;
	timer swatch;
	
	sizes.push_back(8);
	sizes.push_back(10);
	sizes.push_back(12);
	sizes.push_back(14);
	
	for (vector<unsigned>::iterator i=sizes.begin(); i!=sizes.end(); ++i) {
		int count = 1;
		swatch.start();
		result += vandermonde_det(*i);
		// correct for very small times:
		while (swatch.read()<0.02) {
			vandermonde_det(*i);
			++count;
		}
		times.push_back(swatch.read()/count);
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
	return time_vandermonde();
}
