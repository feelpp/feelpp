/** @file time_lw_Qprime.cpp
 *
 *  Test Q' from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
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

static const bool do_test = true;  // set to true in order to run this beast

static unsigned test()
{
	// same matrix as in test P':
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
	matrix m2(m);
	ex a;
	for (unsigned r=0; r<=n*n; ++r) {
		a = m2(r,0);
		for (unsigned c=0; c<n*n; ++c)
			m2.set(r,c,m2(r,c+1));
		m2.set(r,100,a);
	}
	for (unsigned r=0; r<=n*n; ++r)
		for (unsigned c=0; c<=n*n; ++c)
			if (!m(r,c).is_zero())
				m2.set(r,c,m(r,c));
	
	symbol lambda("lambda");
	ex cp = m2.charpoly(lambda);
	
	if (cp.coeff(lambda,0) != numeric("140816284877507872414776")) {
		clog << "characteristic polynomial miscalculated as " << cp << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_Qprime()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;
	
	cout << "timing Lewis-Wester test Q' (charpoly(P'))" << flush;
	
	if (do_test) {
		rolex.start();
		// correct for very small times:
		do {
			result = test();
			++count;
		} while ((time=rolex.read())<0.1 && !result);
		cout << '.' << flush;
		cout << time/count << 's' << endl;
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
	return time_lw_Qprime();
}
