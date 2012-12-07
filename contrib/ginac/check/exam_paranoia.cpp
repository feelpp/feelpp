/** @file exam_paranoia.cpp
 *
 *  This set of tests checks for some of GiNaC's oopses which showed up during
 *  development.  Things were evaluated wrongly and so.  Such a sick behaviour
 *  shouldn't occur any more.  But we are paranoic and we want to exclude these
 *  these oopses for good, so we run those stupid tests... */

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

// The very first pair of historic problems had its roots in power.cpp and was
// finally resolved on April 27th 1999. (Fixing the first on April 23rd
// actually introduced the second.)
static unsigned exam_paranoia1()
{
	unsigned result = 0;
	symbol x("x"), y("y"), z("z");
	ex e, f, g;

	e = x * y * z;
	f = y * z;
	g = e / f;

	// In the first one expand did not do any job at all:
	if (!g.expand().is_equal(x)) {
		clog << "e = x*y*z; f = y*z; expand(e/f) erroneously returned "
		     << g.expand() << endl;
		++result;
	}

	// This one somehow used to return 0:
	e = pow(x + 1, -1);
	if (!e.expand().is_equal(e)) {
		clog << "expand(pow(x + 1, -1)) erroneously returned "
		     << e.expand() << endl;
		++result;
	}

	return result;
}

// And here the second oops which showed up until May 17th 1999.  It had to do
// with lexicographic canonicalization and thus showed up only if the variables
// had the names as given here:
static unsigned exam_paranoia2()
{
	unsigned result = 0;
	symbol x("x"), y("y"), z("z");
	ex e, f, g;

	e = x + z*x;
	f = e*y;
	g = f - e*y;

	// After .expand(), g should be zero:
	if (!g.expand().is_zero()) {
		clog << "e = (x + z*x); f = e*y; expand(f - e*y) erroneously returned "
		     << g.expand() << endl;
		++result;
	}
	// After .eval(), g should be zero:
	if (!g.eval().is_zero()) {
		clog << "e = (x + z*x); f = e*y; eval(f - e*y) erroneously returned "
		     << g.eval() << endl;
		++result;
	}
	// This actually worked already back in April 1999.
	// But we are *very* paranoic!
	if (!g.expand().eval().is_zero()) {
		clog << "e = (x + z*x); f = e*y; eval(expand(f - e*y)) erroneously returned "
		     << g.expand().eval() << endl;
		++result;
	}

	return result;
}

// The third bug was introduced on May 18th 1999, discovered on May 19 and
// fixed that same day.  It worked when x was substituted by 1 but not with
// other numbers:
static unsigned exam_paranoia3()
{
	unsigned result = 0;
	symbol x("x"), y("y");
	ex e, f;

	e = x*y - y;
	f = e.subs(x == 2);

	if (!f.is_equal(y)) {
		clog << "e = x*y - y; f = e.subs(x == 2) erroneously returned "
		     << f << endl;
		++result;
	}
	if (!f.eval().is_equal(y)) {
		clog << "e = x*y - y; eval(e.subs(x == 2)) erroneously returned "
		     << f.eval() << endl;
		++result;
	}
	if (!f.expand().is_equal(y)) {
		clog << "e = x*y - y; expand(e.subs(x == 2)) erroneously returned "
		     << f.expand() << endl;
		++result;
	}

	return result;
}

// The fourth bug was also discovered on May 19th 1999 and fixed immediately:
static unsigned exam_paranoia4()
{
	unsigned result = 0;
	symbol x("x");
	ex e, f, g;

	e = pow(x, 2) + x + 1;
	f = pow(x, 2) + x + 1;
	g = e - f;

	if (!g.is_zero()) {
		clog << "e = pow(x,2) + x + 1; f = pow(x,2) + x + 1; g = e-f; g erroneously returned "
		     << g << endl;
		++result;
	}
	if (!g.is_zero()) {
		clog << "e = pow(x,2) + x + 1; f = pow(x,2) + x + 1; g = e-f; g.eval() erroneously returned "
		     << g.eval() << endl;
		++result;
	}

	return result;
}

// The fifth oops was discovered on May 20th 1999 and fixed a day later:
static unsigned exam_paranoia5()
{
	unsigned result = 0;
	symbol x("x"), y("y");

	ex e, f;
	e = pow(x*y + 1, 2);
	f = pow(x, 2) * pow(y, 2) + 2*x*y + 1;

	if (!(e-f).expand().is_zero()) {
		clog << "e = pow(x*y+1,2); f = pow(x,2)*pow(y,2) + 2*x*y + 1; (e-f).expand() erroneously returned "
		     << (e-f).expand() << endl;
		++result;
	}

	return result;
}

// This one was discovered on Jun 1st 1999 and fixed the same day:
static unsigned exam_paranoia6()
{
	unsigned result = 0;
	symbol x("x");

	ex e, f;
	e = pow(x, -5);
	f = e.denom();

	if (!f.is_equal(pow(x, 5))) {
		clog << "e = pow(x, -5); f = e.denom(); f was " << f << " (should be x^5)" << endl;
		++result;
	}
	return result;
}

// This one was introduced on June 1st 1999 by some aggressive manual
// optimization. Discovered and fixed on June 2nd.
static unsigned exam_paranoia7()
{
	unsigned result = 0;
	symbol x("x"), y("y");

	ex e = y + y*x + 2;
	ex f = expand(pow(e, 2) - (e*y*(x + 1)));

	if (f.nops() > 3) {
		clog << "e=y+y*x+2; f=expand(pow(e,2)-(e*y*(x+1))) has "
		     << f.nops() << " arguments instead of 3 ( f=="
		     << f << " )" << endl;
		++result;
	}
	return result;
}

// This one was a result of the rewrite of mul::max_coefficient when we
// introduced the overall_coefficient field in expairseq objects on Oct 1st
// 1999. Fixed on Oct 4th.
static unsigned exam_paranoia8()
{
	unsigned result = 0;
	symbol x("x");

	ex e = -x / (x+1);
	ex f;
	
	try {
		f = e.normal();
		if (!f.is_equal(e)) {
			clog << "normal(-x/(x+1)) returns " << f << " instead of -x/(x+1)\n";
			++result;
		}
	} catch (const exception &err) {
		clog << "normal(-x/(x+1) throws " << err.what() << endl;
		++result;
	}
	return result;
}

// This one was a result of a modification to frac_cancel() & Co. to avoid
// expanding the numerator and denominator when bringing them from Q[X] to
// Z[X]. multiply_lcm() forgot to multiply the x-linear term with the LCM of
// the coefficient's denominators (2 in this case).  Introduced on Jan 25th
// 2000 and fixed on Jan 31th.
static unsigned exam_paranoia9()
{
	unsigned result = 0;
	symbol x("x");

	ex e = (exp(-x)-2*x*exp(-x)+pow(x,2)/2*exp(-x))/exp(-x);
	ex f = e.normal();

	if (!f.is_equal(1-2*x+pow(x,2)/2)) {
		clog << "normal(" << e << ") returns " << f << " instead of 1-2*x+1/2*x^2\n";
		++result;
	}
	return result;
}

// I have no idea when this broke.  It has been working long ago, before 0.4.0
// and on Feb 13th 2000 I found out that things like 2^(3/2) throw an exception
// "power::eval(): pow(0,0) is undefined" instead of simplifying to 2*2^(1/2).
// It was fixed that same day.
static unsigned exam_paranoia10()
{
	unsigned result = 0;
	
	ex b = numeric(2);
	ex e = numeric(3,2);
	ex r;
	
	try {
		r = pow(b,e).eval();
		if (!(r-2*sqrt(ex(2))).is_zero()) {
			clog << "2^(3/2) erroneously returned " << r << " instead of 2*sqrt(2)" << endl;
			++result;
		}
	} catch (const exception &err) {
		clog << "2^(3/2) throws " << err.what() << endl;
		++result;
	}
	return result;
}

// After the rewriting of basic::normal() & Co. to return {num, den} lists,
// add::normal() forgot to multiply the denominator of the overall_coeff of
// its expanded and normalized children with the denominator of the expanded
// child (did you get this? Well, never mind...). Fixed on Feb 21th 2000.
static unsigned exam_paranoia11()
{
	unsigned result = 0;
	symbol x("x");

	ex e = ((-5-2*x)-((2-5*x)/(-2+x))*(3+2*x))/(5-4*x);
	ex f = e.normal();
	ex d = normal((4+10*x+8*pow(x,2))/(x-2)/(5-4*x));

	if (!(f - d).expand().is_zero()) {
		clog << "normal(" << e << ") returns " << f << " instead of " << d << endl;
		++result;
	}
	return result;
}

// This one returned 0 because add::normal() incorrectly assumed that if the
// common denominator is 1, all the denominators would be 1 (they can in fact
// be +/-1). Fixed on Aug 2nd 2000.
static unsigned exam_paranoia12()
{
	unsigned result = 0;
	symbol x("x");
	
	ex e = 2-2*(1+x)/(-1-x);
	ex f = e.normal();
	ex d = 4;
	
	if (!(f - d).expand().is_zero()) {
		clog << "normal(" << e << ") returns " << f
		     << " instead of " << d << endl;
		++result;
	}
	return result;
}

// This one caused a division by 0 because heur_gcd() didn't check its
// input polynomials against 0. Fixed on Aug 4th 2000.
static unsigned exam_paranoia13()
{
	unsigned result = 0;
	symbol a("a"), b("b"), c("c");
	
	ex e = (b*a-c*a)/(4-a);
	ex d = (c*a-b*a)/(a-4);
	
	try {
		ex f = e.normal();	
		if (!(f - d).expand().is_zero()) {
			clog << "normal(" << e << ") returns " << f
			     << " instead of " << d << endl;
			++result;
		}
	} catch (const exception &err) {
		clog << "normal(" << e << ") throws " << err.what() << endl;
		++result;
	}
	return result;
}

// A bug introduced on July 19, 2001. quo() and rem() would sometimes call
// vector::reserve() with a negative argument. Fixed on Dec 20, 2001.
static unsigned exam_paranoia14()
{
	unsigned result = 0;
	symbol x("x");

	ex q = quo(1, pow(x, 3), x);
	if (!q.is_zero()) {
		clog << "quo(1,x^3,x) erroneously returned " << q << " instead of 0\n";
		++result;
	}

	return result;
}

// Under certain conditions, power::expand_add_2() could produce non-canonical
// numeric expairs. Fixed on Oct 24, 2002.
static unsigned exam_paranoia15()
{
	unsigned result = 0;

	ex q = (pow(pow(2, numeric(1, 2))*2+1, 2)).expand();
	// this used to produce "1+4*sqrt(2)+4*2" which would never evaluate
	// to "9+4*sqrt(2)"

	if (!(q-9-4*pow(2, numeric(1, 2))).is_zero()) {
		clog << "expand((sqrt(2)*2+1)^2) erroneously returned " << q << " instead of 9-4*sqrt(2)\n";
		++result;
	}

	return result;
}

// Expanding products containing powers of sums could return results that
// were not fully expanded. Fixed on Dec 10, 2003.
static unsigned exam_paranoia16()
{
	unsigned result = 0;
	symbol a("a"), b("b"), c("c"), d("d"), e("e");
	ex e1, e2, e3;

	e1 = pow(1+a*sqrt(b+c), 2);
	e2 = e1.expand();

	if (e2.has(pow(a, 2)*(b+c))) {
		clog << "expand(" << e1 << ") didn't fully expand\n";
		++result;
	}

	e1 = (d*sqrt(a+b)+a*sqrt(c+d))*(b*sqrt(a+b)+a*sqrt(c+d));
	e2 = e1.expand();

	if (e2.has(pow(a, 2)*(c+d))) {
		clog << "expand(" << e1 << ") didn't fully expand\n";
		++result;
	}

	e1 = (a+sqrt(b+c))*sqrt(b+c)*(d+sqrt(b+c));
	e2 = e1.expand();

	if (e2.has(a*(b+c))) {
		clog << "expand(" << e1 << ") didn't fully expand\n";
		++result;
	}

	e1 = pow(sqrt(a+b)+sqrt(c+d), 3);
	e2 = e1.expand();

	if (e2.has(3*(a+b)*sqrt(c+d)) || e2.has(3*(c+d)*sqrt(a+b))) {
		clog << "expand(" << e1 << ") didn't fully expand\n";
		++result;
	}

	e1 = a*(b+c*(d+e));
	e2 = e1.expand();

	if (e2.has(c*(d+e))) {
		clog << "expand(" << e1 << ") didn't fully expand\n";
		++result;
	}

	e1 = 2*pow(1+a, 2)/a;
	e2 = e1.expand();

	if (e2.has(pow(a, 2))) {
		clog << "expand(" << e1 << ") didn't fully expand\n";
		++result;
	}

	e1 = a*(a+b);
	e2 = pow(pow(e1, -1), -1);

	if (e2.has(a*b)) {
		clog << "double reciprocal expanded where it should not\n";
		++result;
	}

	return result;
}

// Bug in reposition_dummy_indices() could result in correct expression
// turned into one with inconsistent indices. Fixed on Aug 29, 2006
static unsigned exam_paranoia17()
{
	varidx mu1(symbol("mu1"), 4);
	varidx mu2(symbol("mu2"), 4);
	varidx mu3(symbol("mu3"), 4);
	varidx mu4(symbol("mu4"), 4);
	varidx mu5(symbol("mu5"), 4);
	varidx mu6(symbol("mu6"), 4);

	exvector ev2;
	ev2.push_back(mu3.toggle_variance());
	ev2.push_back(mu6);
	ev2.push_back(mu5.toggle_variance());
	ev2.push_back(mu6.toggle_variance());
	ev2.push_back(mu5);
	ev2.push_back(mu3); 
	// notice: all indices are contracted ...

	ex test_cycl = indexed(symbol("A"), sy_cycl(), ev2);
	test_cycl = test_cycl.simplify_indexed();
	// ... so there should be zero free indices in the end.
	return test_cycl.get_free_indices().size();
}

// Bug in add::eval() could result in numeric terms not being collected into
// the overall coefficient. Fixed on Sep 22, 2010
static unsigned exam_paranoia18()
{
	ex sqrt2 = sqrt(ex(2));
	ex e = 1+2*(sqrt2+1)*(sqrt2-1);
	if ( e.real_part() != 3 ) {
		clog << "real_part(1+2*(sqrt(2)+1)*(sqrt(2)-1)) failed to evaluate to 3\n";
		return 1;
	}
	return 0;
}

// Bug in mul::conjugate when factors are evaluated at branch cuts, reported as
// Sage bug #10964.
static unsigned exam_paranoia19()
{
	symbol a("a");
	ex e = conjugate(a*sqrt(ex(-2))*sqrt(ex(-3)));
	ex c = a*conjugate(sqrt(ex(-2)))*conjugate(sqrt(ex(-3)));
	if (!subs(e-c, a==42).is_zero()) {
		clog << "subs(a*conjugate(sqrt(-2))*conjugate(sqrt(-3))-conjugate(a*sqrt(-2)*sqrt(-3)),a==42) failed to evaluate to 0\n";
		return 1;
	}
	return 0;
}

// Bug in expairseq::is_polynomial (fixed 2011-05-20).
static unsigned exam_paranoia20()
{
	unsigned result = 0;
	symbol x("x");
	ex e1 = sqrt(x*x+1)*sqrt(x+1);
	if (e1.is_polynomial(x)) {
		clog << "sqrt(x*x+1)*sqrt(x+1) is wrongly reported to be a polynomial in x\n";
		++result;
	}
	ex e2 = sqrt(Pi)*x;
	if (!e2.is_polynomial(x)) {
		clog << "sqrt(Pi)*x is wrongly reported to be no polynomial in x\n";
		++result;
	}
	return result;
}

unsigned exam_paranoia()
{
	unsigned result = 0;
	
	cout << "examining several historic failures just out of paranoia" << flush;
	
	result += exam_paranoia1();  cout << '.' << flush;
	result += exam_paranoia2();  cout << '.' << flush;
	result += exam_paranoia3();  cout << '.' << flush;
	result += exam_paranoia4();  cout << '.' << flush;
	result += exam_paranoia5();  cout << '.' << flush;
	result += exam_paranoia6();  cout << '.' << flush;
	result += exam_paranoia7();  cout << '.' << flush;
	result += exam_paranoia8();  cout << '.' << flush;
	result += exam_paranoia9();  cout << '.' << flush;
	result += exam_paranoia10();  cout << '.' << flush;
	result += exam_paranoia11();  cout << '.' << flush;
	result += exam_paranoia12();  cout << '.' << flush;
	result += exam_paranoia13();  cout << '.' << flush;
	result += exam_paranoia14();  cout << '.' << flush;
	result += exam_paranoia15();  cout << '.' << flush;
	result += exam_paranoia16();  cout << '.' << flush;
	result += exam_paranoia17();  cout << '.' << flush;
	result += exam_paranoia18();  cout << '.' << flush;
	result += exam_paranoia19();  cout << '.' << flush;
	result += exam_paranoia20();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_paranoia();
}
