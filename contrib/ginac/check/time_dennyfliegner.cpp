/** @file time_dennyfliegner.cpp
 *
 *  The first test routine implements Denny Fliegner's quick consistency check:
 *  1)  e = (a0 + a1 + a2 + a3 + ...)^2, in expanded form
 *  2)  substitute a0 by (-a2 - a3 - ...) in e
 *  3)  expand e
 *  after which e should be just a1^2. */

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
#include <sstream>
#include <vector>
using namespace std;

static unsigned expand_subs(unsigned size)
{
	unsigned result = 0;
	// create a vector of size symbols named "a0", "a1", ...
	vector<symbol> a;
	ex e;
	for (unsigned i=0; i<size; ++i) {
		ostringstream buf;
		buf << "a" << i << ends;
		a.push_back(symbol(buf.str()));
		e += a[i];
	}
	ex aux;
	
	// prepare aux so it will swallow anything but a1^2:
	aux = -e + a[0] + a[1];
	e = pow(e,2).expand().subs(a[0]==aux, subs_options::no_pattern).expand();
	
	if (e != pow(a[1],2)) {
		clog << "Denny Fliegner's quick consistency check erroneously returned "
		     << e << "." << endl;
		++result;
	}
	
	return result;
}

unsigned time_dennyfliegner()
{
	unsigned result = 0;
	
	cout << "timing commutative expansion and substitution" << flush;
	
	vector<unsigned> sizes;
	vector<double> times;
	timer breitling;
	
	sizes.push_back(100);
	sizes.push_back(200);
	sizes.push_back(400);
	sizes.push_back(800);
	
	for (vector<unsigned>::iterator i=sizes.begin(); i!=sizes.end(); ++i) {
		breitling.start();
		result += expand_subs(*i);
		times.push_back(breitling.read());
		cout << '.' << flush;
	}
	
	// print the report:
	cout << endl << "	size:  ";
	for (vector<unsigned>::iterator i=sizes.begin(); i!=sizes.end(); ++i)
		cout << '\t' << *i;
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
	return time_dennyfliegner();
}
