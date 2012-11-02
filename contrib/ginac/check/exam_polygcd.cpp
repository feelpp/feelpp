/** @file exam_polygcd.cpp
 *
 *  Some test with polynomial GCD calculations. See also the checks for
 *  rational function normalization in normalization.cpp. */

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

const int MAX_VARIABLES = 3;

static symbol x("x"), z("z");
static symbol y[MAX_VARIABLES];

// GCD = 1
static unsigned poly_gcd1()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = x;
		ex e2 = pow(x, 2);
		for (int i=0; i<v; i++) {
			e1 += y[i];
			e2 += pow(y[i], 2);
		}

		ex f = (e1 + 1) * (e1 + 2);
		ex g = e2 * (-pow(x, 2) * y[0] * 3 + pow(y[0], 2) - 1);
		ex r = gcd(f, g);
		if (r != 1) {
			clog << "case 1, gcd(" << f << "," << g << ") = " << r << " (should be 1)" << endl;
			return 1;
		}
	}
	return 0;
}

// Linearly dense quartic inputs with quadratic GCDs
static unsigned poly_gcd2()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = x;
		ex e2 = x;
		for (int i=0; i<v; i++) {
			e1 += y[i];
			e2 -= y[i];
		}

		ex d = pow(e1 + 1, 2);
		ex f = d * pow(e2 - 2, 2);
		ex g = d * pow(e1 + 2, 2);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 2, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Sparse GCD and inputs where degrees are proportional to the number of variables
static unsigned poly_gcd3()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = pow(x, v + 1);
		for (int i=0; i<v; i++)
			e1 += pow(y[i], v + 1);

		ex d = e1 + 1;
		ex f = d * (e1 - 2);
		ex g = d * (e1 + 2);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 3, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Variation of case 3; major performance degradation with PRS
static unsigned poly_gcd3p()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = pow(x, v + 1);
		ex e2 = pow(x, v);
		for (int i=0; i<v; i++) {
			e1 += pow(y[i], v + 1);
			e2 += pow(y[i], v);
		}

		ex d = e1 + 1;
		ex f = d * (e1 - 2);
		ex g = d * (e2 + 2);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 3p, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Quadratic non-monic GCD; f and g have other quadratic factors
static unsigned poly_gcd4()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = pow(x, 2) * pow(y[0], 2);
		ex e2 = pow(x, 2) - pow(y[0], 2);
		ex e3 = x * y[0];
		for (int i=1; i<v; i++) {
			e1 += pow(y[i], 2);
			e2 += pow(y[i], 2);
			e3 += y[i];
		}

		ex d = e1 + 1;
		ex f = d * (e2 - 1);
		ex g = d * pow(e3 + 2, 2);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 4, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Completely dense non-monic quadratic inputs with dense non-monic linear GCDs
static unsigned poly_gcd5()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = x + 1;
		ex e2 = x - 2;
		ex e3 = x + 2;
		for (int i=0; i<v; i++) {
			e1 *= y[i] + 1;
			e2 *= y[i] - 2;
			e3 *= y[i] + 2;
		}

		ex d = e1 - 3;
		ex f = d * (e2 + 3);
		ex g = d * (e3 - 3);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 5, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Sparse non-monic quadratic inputs with linear GCDs
static unsigned poly_gcd5p()
{
	for (int v=1; v<=MAX_VARIABLES; v++) {
		ex e1 = x;
		for (int i=0; i<v; i++)
			e1 *= y[i];

		ex d = e1 - 1;
		ex f = d * (e1 + 3);
		ex g = d * (e1 - 3);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 5p, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Trivariate inputs with increasing degrees
static unsigned poly_gcd6()
{
	symbol y("y");

	for (int j=1; j<=MAX_VARIABLES; j++) {
		ex d = pow(x, j) * y * (z - 1);
		ex f = d * (pow(x, j) + pow(y, j + 1) * pow(z, j) + 1);
		ex g = d * (pow(x, j + 1) + pow(y, j) * pow(z, j + 1) - 7);
		ex r = gcd(f.expand(), g.expand());
		if (!(r - d).expand().is_zero()) {
			clog << "case 6, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
			return 1;
		}
	}
	return 0;
}

// Trivariate polynomials whose GCD has common factors with its cofactors
static unsigned poly_gcd7()
{
	symbol y("y");
	ex p = x - y * z + 1;
	ex q = x - y + z * 3;

	for (int j=1; j<=MAX_VARIABLES; j++) {
		for (int k=j+1; k<=4; k++) {
			ex d = pow(p, j) * pow(q, j);
			ex f = pow(p, j) * pow(q, k);
			ex g = pow(p, k) * pow(q, j); 
			ex r = gcd(f, g);
			if (!(r - d).expand().is_zero() && !(r + d).expand().is_zero()) {
				clog << "case 7, gcd(" << f << "," << g << ") = " << r << " (should be " << d << ")" << endl;
				return 1;
			}
		}
	}
	return 0;
}

unsigned exam_polygcd()
{
	unsigned result = 0;
	
	cout << "examining polynomial GCD computation" << flush;
	
	result += poly_gcd1();  cout << '.' << flush;
	result += poly_gcd2();  cout << '.' << flush;
	result += poly_gcd3();  cout << '.' << flush;
	result += poly_gcd3p();	 cout << '.' << flush; // PRS "worst" case
	result += poly_gcd4();  cout << '.' << flush;
	result += poly_gcd5();  cout << '.' << flush;
	result += poly_gcd5p();  cout << '.' << flush;
	result += poly_gcd6();  cout << '.' << flush;
	result += poly_gcd7();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_polygcd();
}
