/** @file exam_factor.cpp
 *
 *  Factorization test suite. */

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

static unsigned check_factor(const ex& e)
{
	ex ee = e.expand();
	ex answer = factor(ee);
	if ( answer.expand() != ee || answer != e ) {
		clog << "factorization of " << e << " == " << ee << " gave wrong result: " << answer << endl;
		return 1;
	}
	return 0;
}

static unsigned exam_factor1()
{
	unsigned result = 0;
	ex e;
	symbol x("x");
	lst syms;
	syms.append(x);
	
	e = ex("1+x-x^3", syms);
	result += check_factor(e);

	e = ex("1+x^6+x", syms);
	result += check_factor(e);

	e = ex("1-x^6+x", syms);
	result += check_factor(e);

	e = ex("(1+x)^3", syms);
	result += check_factor(e);

	e = ex("(x+1)*(x+4)", syms);
	result += check_factor(e);

	e = ex("x^6-3*x^5+x^4-3*x^3-x^2-3*x+1", syms);
	result += check_factor(e);

	e = ex("(-1+x)^3*(1+x)^3*(1+x^2)", syms);
	result += check_factor(e);

	e = ex("-(-168+20*x-x^2)*(30+x)", syms);
	result += check_factor(e);

	e = ex("x^2*(x-3)^2*(x^3-5*x+7)", syms);
	result += check_factor(e);

	e = ex("-6*x^2*(x-3)", syms);
	result += check_factor(e);

	e = ex("x^16+11*x^4+121", syms);
	result += check_factor(e);

	e = ex("x^8-40*x^6+352*x^4-960*x^2+576", syms);
	result += check_factor(e);

	e = ex("x*(2+x^2)*(1+x+x^3+x^2+x^6+x^5+x^4)*(1+x)^2*(1-x+x^2)^2*(-1+x)", syms);
	result += check_factor(e);

	e = ex("(x+4+x^2-x^3+43*x^4)*(x+1-x^2-3*x^3+4*x^4)", syms);
	result += check_factor(e);

	e = ex("-x^2*(x-1)*(1+x^2)", syms);
	result += check_factor(e);

	e = x;
	result += check_factor(e);

	// x^37 + 1
	e = ex("(1+x)*(1+x^2-x^29-x^11-x^25-x^9-x^35+x^20-x^3+x^16-x^15-x-x^13+x^28+x^24-x^33+x^8-x^19+x^36+x^12-x^27+x^10-x^23+x^18+x^14+x^34-x^31+x^32+x^30-x^5+x^26+x^4+x^22-x^21-x^7-x^17+x^6)", syms);
	result += check_factor(e);

	e = ex("(1+4*x)*x^2*(1-4*x+16*x^2)*(3+5*x+92*x^3)", syms);
	result += check_factor(e);

	e = ex("(77+11*x^3+25*x^2+27*x+102*x^4)*(85+57*x^3+92*x^2+29*x+66*x^4)", syms);
	result += check_factor(e);

	return result;
}

static unsigned exam_factor2()
{
	unsigned result = 0;
	ex e;
	symbol x("x"), y("y"), z("z");
	lst syms;
	syms = x, y, z;
	
	e = ex("x+y", syms);
	result += check_factor(e);

	e = ex("(x^2-y+1)*(x+y)", syms);
	result += check_factor(e);

	e = ex("-2*(x+y)*(x-y)", syms);
	result += check_factor(e);

	e = ex("(16+x^2*z^3)*(-17+3*x-5*z)*(2*x+3*z)*(x-y^2-z^3)", syms);
	result += check_factor(e);

	e = ex("(x-y*z)*(x-y^2-z^3)*(x+y+z)", syms);
	result += check_factor(e);
	
	e = ex("-(y^2-x+z^3)*x*(x+y+z)", syms);
	result += check_factor(e);
	
	e = ex("-316*(3*x-4*z)*(2*x+3*z)*(x+y)*(-1+x)", syms);
	result += check_factor(e);
	
	e = ex("(x+x^3+z^2)*(3*x-4*z)", syms);
	result += check_factor(e);
	
	e = ex("250*(-3+x)*(4*z-3*x)*(x^3+z^2+x)*x", syms);
	result += check_factor(e);
	
	e = ex("327*(x+z^2+x^3)*(3*x-4*z)*(-7+5*x-x^3)*(1+x+x^2)", syms);
	result += check_factor(e);
	
	e = ex("x-y^2-z^3", syms);
	result += check_factor(e);
	
	e = ex("-390*(7+3*x^4)*(2+x^2)*(x-z^3-y^2)", syms);
	result += check_factor(e);
	
	e = ex("55*(1+x)^2*(3*x-4*z)*(1+x+x^2)*(x+x^3+z^2)", syms);
	result += check_factor(e);
	
	e = ex("x+y*x-1", syms);
	result += check_factor(e);
	
	e = ex("390*(-1+x^6-x)*(7+3*x^4)*(2+x^2)*(y+x)*(-1+y-x^2)*(1+x^2+x)^2", syms);
	result += check_factor(e);
	
	e = ex("310*(y+x)*(-1+y-x^2)", syms);
	result += check_factor(e);

	return result;
}

static unsigned exam_factor3()
{
	unsigned result = 0;
	ex e;
	symbol k("k"), n("n");
	lst syms;
	syms = k, n;
	
	e = ex("1/2*(-3+3*k-n)*(-2+3*k-n)*(-1+3*k-n)", syms);
	result += check_factor(e);

	e = ex("1/4*(2*k-n)*(-1+2*k-n)", syms);
	result += check_factor(e);

	return result;
}

unsigned exam_factor()
{
	unsigned result = 0;

	cout << "examining polynomial factorization" << flush;

	result += exam_factor1(); cout << '.' << flush;
	result += exam_factor2(); cout << '.' << flush;
	result += exam_factor3(); cout << '.' << flush;

	return result;
}

int main(int argc, char** argv)
{
	return exam_factor();
}
