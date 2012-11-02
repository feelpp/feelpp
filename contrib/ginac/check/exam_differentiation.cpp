/** @file exam_differentiation.cpp
 *
 *  Tests for symbolic differentiation, including various functions. */

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

static unsigned check_diff(const ex &e, const symbol &x,
						   const ex &d, unsigned nth=1)
{
	ex ed = e.diff(x, nth);
	if (!(ed - d).is_zero()) {
		switch (nth) {
		case 0:
			clog << "zeroth ";
			break;
		case 1:
			break;
		case 2:
			clog << "second ";
			break;
		case 3:
			clog << "third ";
			break;
		default:
			clog << nth << "th ";
		}
		clog << "derivative of " << e << " by " << x << " returned "
		     << ed << " instead of " << d << endl;
		clog << "returned:" << endl;
		clog << tree << ed << "instead of\n" << d << dflt;

		return 1;
	}
	return 0;
}

// Simple (expanded) polynomials
static unsigned exam_differentiation1()
{
	unsigned result = 0;
	symbol x("x"), y("y");
	ex e1, e2, e, d;
	
	// construct bivariate polynomial e to be diff'ed:
	e1 = pow(x, -2) * 3 + pow(x, -1) * 5 + 7 + x * 11 + pow(x, 2) * 13;
	e2 = pow(y, -2) * 5 + pow(y, -1) * 7 + 11 + y * 13 + pow(y, 2) * 17;
	e = (e1 * e2).expand();
	
	// d e / dx:
	d = ex("121-55/x^2-66/x^3-30/x^3/y^2-42/x^3/y-78/x^3*y-102/x^3*y^2-25/x^2/y^2-35/x^2/y-65/x^2*y-85/x^2*y^2+77/y+143*y+187*y^2+130*x/y^2+182/y*x+338*x*y+442*x*y^2+55/y^2+286*x",lst(x,y));
	result += check_diff(e, x, d);
	
	// d e / dy:
	d = ex("91-30/x^2/y^3-21/x^2/y^2+39/x^2+102/x^2*y-50/x/y^3-35/x/y^2+65/x+170/x*y-77*x/y^2+143*x+374*x*y-130/y^3*x^2-91/y^2*x^2+169*x^2+442*x^2*y-110/y^3*x-70/y^3+238*y-49/y^2",lst(x,y));
	result += check_diff(e, y, d);
	
	// d^2 e / dx^2:
	d = ex("286+90/x^4/y^2+126/x^4/y+234/x^4*y+306/x^4*y^2+50/x^3/y^2+70/x^3/y+130/x^3*y+170/x^3*y^2+130/y^2+182/y+338*y+442*y^2+198/x^4+110/x^3",lst(x,y));
	result += check_diff(e, x, d, 2);
	
	// d^2 e / dy^2:
	d = ex("238+90/x^2/y^4+42/x^2/y^3+102/x^2+150/x/y^4+70/x/y^3+170/x+330*x/y^4+154*x/y^3+374*x+390*x^2/y^4+182*x^2/y^3+442*x^2+210/y^4+98/y^3",lst(x,y));
	result += check_diff(e, y, d, 2);
	
	return result;
}

// Trigonometric functions
static unsigned exam_differentiation2()
{
	unsigned result = 0;
	symbol x("x"), y("y"), a("a"), b("b");
	ex e1, e2, e, d;
	
	// construct expression e to be diff'ed:
	e1 = y*pow(x, 2) + a*x + b;
	e2 = sin(e1);
	e = b*pow(e2, 2) + y*e2 + a;
	
	d = 2*b*e2*cos(e1)*(2*x*y + a) + y*cos(e1)*(2*x*y + a);
	result += check_diff(e, x, d);
	
	d = 2*b*pow(cos(e1),2)*pow(2*x*y + a, 2) + 4*b*y*e2*cos(e1)
	    - 2*b*pow(e2,2)*pow(2*x*y + a, 2) - y*e2*pow(2*x*y + a, 2)
	    + 2*pow(y,2)*cos(e1);
	result += check_diff(e, x, d, 2);
	
	d = 2*b*e2*cos(e1)*pow(x, 2) + e2 + y*cos(e1)*pow(x, 2);
	result += check_diff(e, y, d);

	d = 2*b*pow(cos(e1),2)*pow(x,4) - 2*b*pow(e2,2)*pow(x,4)
	    + 2*cos(e1)*pow(x,2) - y*e2*pow(x,4);
	result += check_diff(e, y, d, 2);
	
	// construct expression e to be diff'ed:
	e2 = cos(e1);
	e = b*pow(e2, 2) + y*e2 + a;
	
	d = -2*b*e2*sin(e1)*(2*x*y + a) - y*sin(e1)*(2*x*y + a);
	result += check_diff(e, x, d);
	
	d = 2*b*pow(sin(e1),2)*pow(2*y*x + a,2) - 4*b*e2*sin(e1)*y 
	    - 2*b*pow(e2,2)*pow(2*y*x + a,2) - y*e2*pow(2*y*x + a,2)
	    - 2*pow(y,2)*sin(e1);
	result += check_diff(e, x, d, 2);
	
	d = -2*b*e2*sin(e1)*pow(x,2) + e2 - y*sin(e1)*pow(x, 2);
	result += check_diff(e, y, d);
	
	d = -2*b*pow(e2,2)*pow(x,4) + 2*b*pow(sin(e1),2)*pow(x,4)
	    - 2*sin(e1)*pow(x,2) - y*e2*pow(x,4);
	result += check_diff(e, y, d, 2);

	return result;
}
	
// exp function
static unsigned exam_differentiation3()
{
	unsigned result = 0;
	symbol x("x"), y("y"), a("a"), b("b");
	ex e1, e2, e, d;

	// construct expression e to be diff'ed:
	e1 = y*pow(x, 2) + a*x + b;
	e2 = exp(e1);
	e = b*pow(e2, 2) + y*e2 + a;
	
	d = 2*b*pow(e2, 2)*(2*x*y + a) + y*e2*(2*x*y + a);
	result += check_diff(e, x, d);
	
	d = 4*b*pow(e2,2)*pow(2*y*x + a,2) + 4*b*pow(e2,2)*y
	    + 2*pow(y,2)*e2 + y*e2*pow(2*y*x + a,2);
	result += check_diff(e, x, d, 2);
	
	d = 2*b*pow(e2,2)*pow(x,2) + e2 + y*e2*pow(x,2);
	result += check_diff(e, y, d);
	
	d = 4*b*pow(e2,2)*pow(x,4) + 2*e2*pow(x,2) + y*e2*pow(x,4);
	result += check_diff(e, y, d, 2);

	return result;
}

// log functions
static unsigned exam_differentiation4()
{
	unsigned result = 0;
	symbol x("x"), y("y"), a("a"), b("b");
	ex e1, e2, e, d;
	
	// construct expression e to be diff'ed:
	e1 = y*pow(x, 2) + a*x + b;
	e2 = log(e1);
	e = b*pow(e2, 2) + y*e2 + a;
	
	d = 2*b*e2*(2*x*y + a)/e1 + y*(2*x*y + a)/e1;
	result += check_diff(e, x, d);
	
	d = 2*b*pow((2*x*y + a),2)*pow(e1,-2) + 4*b*y*e2/e1
	    - 2*b*e2*pow(2*x*y + a,2)*pow(e1,-2) + 2*pow(y,2)/e1
	    - y*pow(2*x*y + a,2)*pow(e1,-2);
	result += check_diff(e, x, d, 2);
	
	d = 2*b*e2*pow(x,2)/e1 + e2 + y*pow(x,2)/e1;
	result += check_diff(e, y, d);
	
	d = 2*b*pow(x,4)*pow(e1,-2) - 2*b*e2*pow(e1,-2)*pow(x,4)
	    + 2*pow(x,2)/e1 - y*pow(x,4)*pow(e1,-2);
	result += check_diff(e, y, d, 2);

	return result;
}

// Functions with two variables
static unsigned exam_differentiation5()
{
	unsigned result = 0;
	symbol x("x"), y("y"), a("a"), b("b");
	ex e1, e2, e, d;
	
	// test atan2
	e1 = y*pow(x, 2) + a*x + b;
	e2 = x*pow(y, 2) + b*y + a;
	e = atan2(e1,e2);
	
	d = pow(y,2)*pow(pow(b+y*pow(x,2)+x*a,2)+pow(y*b+pow(y,2)*x+a,2),-1)*
	    (-b-y*pow(x,2)-x*a)
	   +pow(pow(b+y*pow(x,2)+x*a,2)+pow(y*b+pow(y,2)*x+a,2),-1)*
	    (y*b+pow(y,2)*x+a)*(2*y*x+a);
	result += check_diff(e, x, d);
	
	return result;
}

// Series
static unsigned exam_differentiation6()
{
	symbol x("x");
	ex e, d, ed;
	
	e = sin(x).series(x==0, 8);
	d = cos(x).series(x==0, 7);
	ed = e.diff(x);
	ed = series_to_poly(ed);
	d = series_to_poly(d);
	
	if (!(ed - d).is_zero()) {
		clog << "derivative of " << e << " by " << x << " returned "
		     << ed << " instead of " << d << ")" << endl;
		return 1;
	}
	return 0;
}

// Hashing can help a lot, if differentiation is done cleverly
static unsigned exam_differentiation7()
{
	symbol x("x");
	ex P = x + pow(x,3);
	ex e = (P.diff(x) / P).diff(x, 2);
	ex d = 6/P - 18*x/pow(P,2) - 54*pow(x,3)/pow(P,2) + 2/pow(P,3)
	    +18*pow(x,2)/pow(P,3) + 54*pow(x,4)/pow(P,3) + 54*pow(x,6)/pow(P,3);
	
	if (!(e-d).expand().is_zero()) {
		clog << "expanded second derivative of " << (P.diff(x) / P) << " by " << x
		     << " returned " << e.expand() << " instead of " << d << endl;
		return 1;
	}
	if (e.nops() > 3) {
		clog << "second derivative of " << (P.diff(x) / P) << " by " << x
		     << " has " << e.nops() << " operands.  "
		     << "The result is still correct but not optimal: 3 are enough!  "
		     << "(Hint: maybe the product rule for objects of class mul should be more careful about assembling the result?)" << endl;
		return 1;
	}
	return 0;
}

unsigned exam_differentiation()
{
	unsigned result = 0;
	
	cout << "examining symbolic differentiation" << flush;
	
	result += exam_differentiation1();  cout << '.' << flush;
	result += exam_differentiation2();  cout << '.' << flush;
	result += exam_differentiation3();  cout << '.' << flush;
	result += exam_differentiation4();  cout << '.' << flush;
	result += exam_differentiation5();  cout << '.' << flush;
	result += exam_differentiation6();  cout << '.' << flush;
	result += exam_differentiation7();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_differentiation();
}
