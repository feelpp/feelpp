/** @File exam_pseries.cpp
 *
 *  Series expansion test (Laurent and Taylor series). */

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

static symbol x("x");

static unsigned check_series(const ex &e, const ex &point, const ex &d, int order = 8)
{
	ex es = e.series(x==point, order);
	ex ep = ex_to<pseries>(es).convert_to_poly();
	if (!(ep - d).expand().is_zero()) {
		clog << "series expansion of " << e << " at " << point
		     << " erroneously returned " << ep << " (instead of " << d
		     << ")" << endl;
		clog << tree << (ep-d) << dflt;
		return 1;
	}
	return 0;
}

// Series expansion
static unsigned exam_series1()
{
	using GiNaC::log;

	symbol a("a");
	symbol b("b");
	unsigned result = 0;
	ex e, d;
	
	e = pow(a+b, x);
	d = 1 + Order(pow(x, 1));
	result += check_series(e, 0, d, 1);

	e = sin(x);
	d = x - pow(x, 3) / 6 + pow(x, 5) / 120 - pow(x, 7) / 5040 + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = cos(x);
	d = 1 - pow(x, 2) / 2 + pow(x, 4) / 24 - pow(x, 6) / 720 + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = exp(x);
	d = 1 + x + pow(x, 2) / 2 + pow(x, 3) / 6 + pow(x, 4) / 24 + pow(x, 5) / 120 + pow(x, 6) / 720 + pow(x, 7) / 5040 + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = pow(1 - x, -1);
	d = 1 + x + pow(x, 2) + pow(x, 3) + pow(x, 4) + pow(x, 5) + pow(x, 6) + pow(x, 7) + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = x + pow(x, -1);
	d = x + pow(x, -1);
	result += check_series(e, 0, d);
	
	e = x + pow(x, -1);
	d = 2 + pow(x-1, 2) - pow(x-1, 3) + pow(x-1, 4) - pow(x-1, 5) + pow(x-1, 6) - pow(x-1, 7) + Order(pow(x-1, 8));
	result += check_series(e, 1, d);
	
	e = pow(x + pow(x, 3), -1);
	d = pow(x, -1) - x + pow(x, 3) - pow(x, 5) + pow(x, 7) + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = pow(pow(x, 2) + pow(x, 4), -1);
	d = pow(x, -2) - 1 + pow(x, 2) - pow(x, 4) + pow(x, 6) + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = pow(sin(x), -2);
	d = pow(x, -2) + numeric(1,3) + pow(x, 2) / 15 + pow(x, 4) * 2/189 + pow(x, 6) / 675  + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = sin(x) / cos(x);
	d = x + pow(x, 3) / 3 + pow(x, 5) * 2/15 + pow(x, 7) * 17/315 + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = cos(x) / sin(x);
	d = pow(x, -1) - x / 3 - pow(x, 3) / 45 - pow(x, 5) * 2/945 - pow(x, 7) / 4725 + Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	e = pow(numeric(2), x);
	ex t = log(2) * x;
	d = 1 + t + pow(t, 2) / 2 + pow(t, 3) / 6 + pow(t, 4) / 24 + pow(t, 5) / 120 + pow(t, 6) / 720 + pow(t, 7) / 5040 + Order(pow(x, 8));
	result += check_series(e, 0, d.expand());
	
	e = pow(Pi, x);
	t = log(Pi) * x;
	d = 1 + t + pow(t, 2) / 2 + pow(t, 3) / 6 + pow(t, 4) / 24 + pow(t, 5) / 120 + pow(t, 6) / 720 + pow(t, 7) / 5040 + Order(pow(x, 8));
	result += check_series(e, 0, d.expand());
	
	e = log(x);
	d = e;
	result += check_series(e, 0, d, 1);
	result += check_series(e, 0, d, 2);
	
	e = pow(x, 8) * pow(pow(x,3)+ pow(x + pow(x,3), 2), -2);
	d = pow(x, 4) - 2*pow(x, 5) + Order(pow(x, 6));
	result += check_series(e, 0, d, 6);
	
	e = cos(x) * pow(sin(x)*(pow(x, 5) + 4 * pow(x, 2)), -3);
	d = pow(x, -9) / 64 - 3 * pow(x, -6) / 256 - pow(x, -5) / 960 + 535 * pow(x, -3) / 96768
	    + pow(x, -2) / 1280 - pow(x, -1) / 14400 - numeric(283, 129024) - 2143 * x / 5322240
	    + Order(pow(x, 2));
	result += check_series(e, 0, d, 2);
	
	e = sqrt(1+x*x) * sqrt(1+2*x*x);
	d = 1 + Order(pow(x, 2));
	result += check_series(e, 0, d, 2);

	e = pow(x, 4) * sin(a) + pow(x, 2);
	d = pow(x, 2) + Order(pow(x, 3));
	result += check_series(e, 0, d, 3);

	e = log(a*x + b*x*x*log(x));
	d = log(a*x) + b/a*log(x)*x - pow(b/a, 2)/2*pow(log(x)*x, 2) + Order(pow(x, 3));
	result += check_series(e, 0, d, 3);

	e = pow((x+a), b);
	d = pow(a, b) + (pow(a, b)*b/a)*x + (pow(a, b)*b*b/a/a/2 - pow(a, b)*b/a/a/2)*pow(x, 2) + Order(pow(x, 3));
	result += check_series(e, 0, d, 3);

	return result;
}

// Series addition
static unsigned exam_series2()
{
	unsigned result = 0;
	ex e, d;
	
	e = pow(sin(x), -1).series(x==0, 8) + pow(sin(-x), -1).series(x==0, 12);
	d = Order(pow(x, 8));
	result += check_series(e, 0, d);
	
	return result;
}

// Series multiplication
static unsigned exam_series3()
{
	unsigned result = 0;
	ex e, d;
	
	e = sin(x).series(x==0, 8) * pow(sin(x), -1).series(x==0, 12);
	d = 1 + Order(pow(x, 7));
	result += check_series(e, 0, d);
	
	return result;
}

// Series exponentiation
static unsigned exam_series4()
{
	using GiNaC::tgamma;
	unsigned result = 0;
	ex e, d;
	
	e = pow((2*cos(x)).series(x==0, 5), 2).series(x==0, 5);
	d = 4 - 4*pow(x, 2) + 4*pow(x, 4)/3 + Order(pow(x, 5));
	result += check_series(e, 0, d);
	
	e = pow(tgamma(x), 2).series(x==0, 2);
	d = pow(x,-2) - 2*Euler/x + (pow(Pi,2)/6+2*pow(Euler,2)) 
		+ x*(-4*pow(Euler, 3)/3 -pow(Pi,2)*Euler/3 - 2*zeta(3)/3) + Order(pow(x, 2));
	result += check_series(e, 0, d);
	
	return result;
}

// Order term handling
static unsigned exam_series5()
{
	unsigned result = 0;
	ex e, d;

	e = 1 + x + pow(x, 2) + pow(x, 3);
	d = Order(1);
	result += check_series(e, 0, d, 0);
	d = 1 + Order(x);
	result += check_series(e, 0, d, 1);
	d = 1 + x + Order(pow(x, 2));
	result += check_series(e, 0, d, 2);
	d = 1 + x + pow(x, 2) + Order(pow(x, 3));
	result += check_series(e, 0, d, 3);
	d = 1 + x + pow(x, 2) + pow(x, 3);
	result += check_series(e, 0, d, 4);
	return result;
}

// Series expansion of tgamma(-1)
static unsigned exam_series6()
{
	using GiNaC::tgamma;
	ex e = tgamma(2*x);
	ex d = pow(x+1,-1)*numeric(1,4) +
	       pow(x+1,0)*(numeric(3,4) -
	                   numeric(1,2)*Euler) +
	       pow(x+1,1)*(numeric(7,4) -
	                   numeric(3,2)*Euler +
	                   numeric(1,2)*pow(Euler,2) +
	                   numeric(1,12)*pow(Pi,2)) +
	       pow(x+1,2)*(numeric(15,4) -
	                   numeric(7,2)*Euler -
	                   numeric(1,3)*pow(Euler,3) +
	                   numeric(1,4)*pow(Pi,2) +
	                   numeric(3,2)*pow(Euler,2) -
	                   numeric(1,6)*pow(Pi,2)*Euler -
	                   numeric(2,3)*zeta(3)) +
	       pow(x+1,3)*(numeric(31,4) - pow(Euler,3) -
	                   numeric(15,2)*Euler +
	                   numeric(1,6)*pow(Euler,4) +
	                   numeric(7,2)*pow(Euler,2) +
	                   numeric(7,12)*pow(Pi,2) -
	                   numeric(1,2)*pow(Pi,2)*Euler -
	                   numeric(2)*zeta(3) +
	                   numeric(1,6)*pow(Euler,2)*pow(Pi,2) +
	                   numeric(1,40)*pow(Pi,4) +
	                   numeric(4,3)*zeta(3)*Euler) +
	       Order(pow(x+1,4));
	return check_series(e, -1, d, 4);
}
	
// Series expansion of tan(x==Pi/2)
static unsigned exam_series7()
{
	ex e = tan(x*Pi/2);
	ex d = pow(x-1,-1)/Pi*(-2) + pow(x-1,1)*Pi/6 + pow(x-1,3)*pow(Pi,3)/360
	      +pow(x-1,5)*pow(Pi,5)/15120 + pow(x-1,7)*pow(Pi,7)/604800
	      +Order(pow(x-1,9));
	return check_series(e,1,d,9);
}

// Series expansion of log(sin(x==0))
static unsigned exam_series8()
{
	ex e = log(sin(x));
	ex d = log(x) - pow(x,2)/6 - pow(x,4)/180 - pow(x,6)/2835 - pow(x,8)/37800 + Order(pow(x,9));
	return check_series(e,0,d,9);
}

// Series expansion of Li2(sin(x==0))
static unsigned exam_series9()
{
	ex e = Li2(sin(x));
	ex d = x + pow(x,2)/4 - pow(x,3)/18 - pow(x,4)/48
	       - 13*pow(x,5)/1800 - pow(x,6)/360 - 23*pow(x,7)/21168
	       + Order(pow(x,8));
	return check_series(e,0,d,8);
}

// Series expansion of Li2((x==2)^2), caring about branch-cut
static unsigned exam_series10()
{
	using GiNaC::log;

	ex e = Li2(pow(x,2));
	ex d = Li2(4) + (-log(3) + I*Pi*csgn(I-I*pow(x,2))) * (x-2)
	       + (numeric(-2,3) + log(3)/4 - I*Pi/4*csgn(I-I*pow(x,2))) * pow(x-2,2)
	       + (numeric(11,27) - log(3)/12 + I*Pi/12*csgn(I-I*pow(x,2))) * pow(x-2,3)
	       + (numeric(-155,648) + log(3)/32 - I*Pi/32*csgn(I-I*pow(x,2))) * pow(x-2,4)
	       + Order(pow(x-2,5));
	return check_series(e,2,d,5);
}

// Series expansion of logarithms around branch points
static unsigned exam_series11()
{
	using GiNaC::log;

	unsigned result = 0;
	ex e, d;
	symbol a("a");
	
	e = log(x);
	d = log(x);
	result += check_series(e,0,d,5);
	
	e = log(3/x);
	d = log(3)-log(x);
	result += check_series(e,0,d,5);
	
	e = log(3*pow(x,2));
	d = log(3)+2*log(x);
	result += check_series(e,0,d,5);
	
	// These ones must not be expanded because it would result in a branch cut
	// running in the wrong direction. (Other systems tend to get this wrong.)
	e = log(-x);
	d = e;
	result += check_series(e,0,d,5);
	
	e = log(I*(x-123));
	d = e;
	result += check_series(e,123,d,5);
	
	e = log(a*x);
	d = e;  // we don't know anything about a!
	result += check_series(e,0,d,5);
	
	e = log((1-x)/x);
	d = log(1-x) - (x-1) + pow(x-1,2)/2 - pow(x-1,3)/3  + pow(x-1,4)/4 + Order(pow(x-1,5));
	result += check_series(e,1,d,5);
	
	return result;
}

// Series expansion of other functions around branch points
static unsigned exam_series12()
{
	using GiNaC::log;
	using GiNaC::atanh;

	unsigned result = 0;
	ex e, d;
	
	// NB: Mma and Maple give different results, but they agree if one
	// takes into account that by assumption |x|<1.
	e = atan(x);
	d = (I*log(2)/2-I*log(1+I*x)/2) + (x-I)/4 + I*pow(x-I,2)/16 + Order(pow(x-I,3));
	result += check_series(e,I,d,3);
	
	// NB: here, at -I, Mathematica disagrees, but it is wrong -- they
	// pick up a complex phase by incorrectly expanding logarithms.
	e = atan(x);
	d = (-I*log(2)/2+I*log(1-I*x)/2) + (x+I)/4 - I*pow(x+I,2)/16 + Order(pow(x+I,3));
	result += check_series(e,-I,d,3);
	
	// This is basically the same as above, the branch point is at +/-1:
	e = atanh(x);
	d = (-log(2)/2+log(x+1)/2) + (x+1)/4 + pow(x+1,2)/16 + Order(pow(x+1,3));
	result += check_series(e,-1,d,3);
	
	return result;
}

// Test of the patch of Stefan Weinzierl that prevents an infinite loop if
// a factor in a product is a complicated way of writing zero.
static unsigned exam_series13()
{
	unsigned result = 0;

	ex e = (new mul(pow(2,x), (1/x*(-(1+x)/(1-x)) + (1+x)/x/(1-x)))
	       )->setflag(status_flags::evaluated);
	ex d = Order(x);
	result += check_series(e,0,d,1);

	return result;
}

// Test if (1+x)^(1/x) can be expanded.
static unsigned exam_series14()
{
	unsigned result = 0;

	ex e = pow(1+x, sin(x)/x);
	ex d = 1 + x - pow(x,3)/6 + Order(pow(x,4));
	try {
		result += check_series(e,0,d,4);
	} catch (const pole_error& err) {
		clog << "series expansion of " << e << " at 0 raised an exception." << endl;
		++result;
	}

	return result;
}

unsigned exam_pseries()
{
	unsigned result = 0;
	
	cout << "examining series expansion" << flush;
	
	result += exam_series1();  cout << '.' << flush;
	result += exam_series2();  cout << '.' << flush;
	result += exam_series3();  cout << '.' << flush;
	result += exam_series4();  cout << '.' << flush;
	result += exam_series5();  cout << '.' << flush;
	result += exam_series6();  cout << '.' << flush;
	result += exam_series7();  cout << '.' << flush;
	result += exam_series8();  cout << '.' << flush;
	result += exam_series9();  cout << '.' << flush;
	result += exam_series10();  cout << '.' << flush;
	result += exam_series11();  cout << '.' << flush;
	result += exam_series12();  cout << '.' << flush;
	result += exam_series13();  cout << '.' << flush;
	result += exam_series14();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_pseries();
}
