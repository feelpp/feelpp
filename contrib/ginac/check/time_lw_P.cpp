/** @file time_lw_P.cpp
 *
 *  Test P from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
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
	// This is a pattern that comes up in graph theory:
	const unsigned n = 10;
	matrix m(n*n+1,n*n+1);
	for (unsigned i=1; i<=n*n; ++i)
		m.set(i-1,i-1,1);
	for (unsigned i=1; i<=n*n; ++i)
		if (!(i%n))
			m.set(i-1,n*n,1);
	for (unsigned i=1; i<=n*n; ++i)
		if (!((i-1)%n))
			m.set(n*n,i-1,n-(i-1)/n);
	for(unsigned i=1; i<=n; ++i)
		for (unsigned j=1; j<=n; ++j)
			if (i-j)
				for (unsigned k=1; k<n; ++k)
					m.set((i-1)*n+k-1,(j-1)*n+k,n+1-j);
	
	ex det = m.determinant();
	
	if (det!=numeric("75810815066186520")) {
		clog << "det of sparse rank 101 matrix erroneously returned " << det << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_P()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;
	
	cout << "timing Lewis-Wester test P (det of sparse rank 101)" << flush;
	
	rolex.start();
	// correct for very small times:
	do {
		result = test();
		++count;
	} while ((time=rolex.read())<0.1 && !result);
	cout << time/count << 's' << endl;
	
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_P();
}
