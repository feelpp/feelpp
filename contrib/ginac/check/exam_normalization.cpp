/** @file exam_normalization.cpp
 *
 *  Rational function normalization test suite. */

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

static symbol w("w"), x("x"), y("y"), z("z");

static unsigned check_normal(const ex &e, const ex &d)
{
	ex en = e.normal();
	if (!en.is_equal(d)) {
		clog << "normal form of " << e << " erroneously returned "
		     << en << " (should be " << d << ")" << endl;
		return 1;
	}
	return 0;
}

static unsigned exam_normal1()
{
	unsigned result = 0;
	ex e, d;
	
	// Expansion
	e = pow(x, 2) - (x+1)*(x-1) - 1;
	d = 0;
	result += check_normal(e, d);
	
	// Expansion inside functions
	e = sin(x*(x+1)-x) + 1;
	d = sin(pow(x, 2)) + 1;
	result += check_normal(e, d);
	
	// Fraction addition
	e = 2/x + y/3;
	d = (x*y + 6) / (x*3);
	result += check_normal(e, d);
	
	e = pow(x, -1) + x/(x+1);
	d = (pow(x, 2)+x+1)/(x*(x+1));
	result += check_normal(e, d);

	return result;
}

static unsigned exam_normal2()
{
	unsigned result = 0;
	ex e, d;
	
	// Fraction cancellation
	e = numeric(1)/2 * z * (2*x + 2*y);
	d = z * (x + y);
	result += check_normal(e, d);
	
	e = numeric(1)/6 * z * (3*x + 3*y) * (2*x + 2*w);
	d = z * (x + y) * (x + w);
	result += check_normal(e, d);
	
	e = (3*x + 3*y) * (w/3 + z/3);
	d = (x + y) * (w + z);
	result += check_normal(e, d);
	
	// Fails stochastically with the new tinfo mechanism, because
	// sometimes the equivalent answer ... / pow(y - x, 2) is calculated.
	// TODO: make check for both cases.
//	e = (pow(x, 2) - pow(y, 2)) / pow(x-y, 3);
//	d = (x + y) / pow(x - y, 2);
//	result += check_normal(e, d);
	
	e = (pow(x, -1) + x) / (pow(x , 2) * 2 + 2);
	d = pow(x * 2, -1);
	result += check_normal(e, d);
	
	// Fails stochastically with the new tinfo mechanism, because
	// sometimes the equivalent answer ... / pow(y - x, 2) is calculated.
	// TODO: make check for both cases.
	// Fraction cancellation with rational coefficients
//	e = (pow(x, 2) - pow(y, 2)) / pow(x/2 - y/2, 3);
//	d = (8 * x + 8 * y) / pow(x - y, 2);
//	result += check_normal(e, d);
	
	// Fraction cancellation with rational coefficients
	e = z/5 * (x/7 + y/10) / (x/14 + y/20);
	d = 2*z/5;
	result += check_normal(e, d);
	
	return result;
}

static unsigned exam_normal3()
{
	unsigned result = 0;
	ex e, d;
	
	// Distribution of powers
	e = pow(x/y, 2);
	d = pow(x, 2) / pow(y, 2);
	result += check_normal(e, d);
	
	// Distribution of powers (integer, distribute) and fraction addition
	e = pow(pow(x, -1) + x, 2);
	d = pow(pow(x, 2) + 1, 2) / pow(x, 2);
	result += check_normal(e, d);
	
	// Distribution of powers (non-integer, don't distribute) and fraction addition
	e = pow(pow(x, -1) + x, numeric(1)/2);
	d = pow((pow(x, 2) + 1) / x, numeric(1)/2);
	result += check_normal(e, d);
	
	return result;
}

static unsigned exam_normal4()
{
	unsigned result = 0;
	ex e, d;
	
	// Replacement of functions with temporary symbols and fraction cancellation
	e = pow(sin(x), 2) - pow(cos(x), 2);
	e /= sin(x) + cos(x);
	d = sin(x) - cos(x);
	result += check_normal(e, d);
	
	// Replacement of non-integer powers with temporary symbols
	e = (pow(numeric(2), numeric(1)/2) * x + x) / x;
	d = pow(numeric(2), numeric(1)/2) + 1;
	result += check_normal(e, d);
	
	// Replacement of complex numbers with temporary symbols
	e = (x + y + x*I + y*I) / (x + y);
	d = 1 + I;
	result += check_normal(e, d);
	
	e = (pow(x, 2) + pow(y, 2)) / (x + y*I);
	d = e;
	result += check_normal(e, d);
	
	// More complex rational function
	e = (pow(x-y*2,4)/pow(pow(x,2)-pow(y,2)*4,2)+1)*(x+y*2)*(y+z)/(pow(x,2)+pow(y,2)*4);
	d = (y*2 + z*2) / (x + y*2);
	result += check_normal(e, d);
	
	return result;
}

/* Test content(), integer_content(), primpart(). */
static unsigned check_content(const ex & e, const ex & x, const ex & ic, const ex & c, const ex & pp)
{
	unsigned result = 0;

	ex r_ic = e.integer_content();
	if (!r_ic.is_equal(ic)) {
		clog << "integer_content(" << e << ") erroneously returned "
		     << r_ic << " instead of " << ic << endl;
		++result;
	}

	ex r_c = e.content(x);
	if (!r_c.is_equal(c)) {
		clog << "content(" << e << ", " << x << ") erroneously returned "
		     << r_c << " instead of " << c << endl;
		++result;
	}

	ex r_pp = e.primpart(x);
	if (!r_pp.is_equal(pp)) {
		clog << "primpart(" << e << ", " << x << ") erroneously returned "
		     << r_pp << " instead of " << pp << endl;
		++result;
	}

	ex r = r_c*r_pp*e.unit(x);
	if (!(r - e).expand().is_zero()) {
		clog << "product of unit, content, and primitive part of " << e << " yielded "
		     << r << " instead of " << e << endl;
		++result;
	}

	return result;
}

static unsigned exam_content()
{
	unsigned result = 0;
	symbol x("x"), y("y");

	result += check_content(ex(-3)/4, x, ex(3)/4, ex(3)/4, 1);
	result += check_content(-x/4, x, ex(1)/4, ex(1)/4, x);
	result += check_content(5*x-15, x, 5, 5, x-3);
	result += check_content(5*x*y-15*y*y, x, 5, 5*y, x-3*y);
	result += check_content(-15*x/2+ex(25)/3, x, ex(5)/6, ex(5)/6, 9*x-10);
	result += check_content(-x*y, x, 1, y, x);

	return result;
}

unsigned exam_normalization()
{
	unsigned result = 0;
	
	cout << "examining rational function normalization" << flush;
	
	result += exam_normal1(); cout << '.' << flush;
	result += exam_normal2(); cout << '.' << flush;
	result += exam_normal3(); cout << '.' << flush;
	result += exam_normal4(); cout << '.' << flush;
	result += exam_content(); cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_normalization();
}
