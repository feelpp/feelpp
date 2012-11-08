/** @file exam_misc.cpp
 *
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
using namespace GiNaC;

#include <iostream>
using namespace std;

#define VECSIZE 30
static unsigned exam_expand_subs()
{
	unsigned result = 0;
	symbol a[VECSIZE];
	ex e, aux;
	
	for (unsigned i=0; i<VECSIZE; ++i)
		e = e + a[i];
	
	// prepare aux so it will swallow anything but a1^2:
	aux = -e + a[0] + a[1];
	e = expand(subs(expand(pow(e, 2)), a[0] == aux));
	
	if (e != pow(a[1],2)) {
		clog << "Denny Fliegner's quick consistency check erroneously returned "
		     << e << "." << endl;
		++result;
	}
	
	return result;
}

/*  A simple modification of Denny Fliegner's three step consistency test:
 *  1)  e = (a0 + a1)^200
 *  2)  expand e
 *  3)  substitute a0 by -a1 in e
 *  after which e should return 0 (without expanding). */
static unsigned exam_expand_subs2()
{
	unsigned result = 0;
	symbol a("a"), b("b");
	ex e, f;
	
	e = pow(a+b,200).expand();
	f = e.subs(a == -b);
	
	if (f != 0) {
		clog << "e = pow(a+b,200).expand(); f = e.subs(a == -b); erroneously returned "
		     << f << " instead of simplifying to 0." << endl;
		++result;
	}
	
	return result;
}

static unsigned exam_expand_power()
{
	unsigned result = 0;
	symbol x("x"), a("a"), b("b");
	ex e;
	
	e = pow(x,pow(a+b,2)-pow(a,2)-pow(b,2)-a*b*2).expand();
	
	if (e != 1) {
		clog << "e = pow(x,pow(a+b,2)-pow(a,2)-pow(b,2)-a*b*2).expand(); erroneously returned "
		     << e << " instead of simplifying to 1." << endl;
		++result;
	}
	
	return result;
}

static unsigned exam_sqrfree()
{
	unsigned result = 0;
	symbol x("x"), y("y");
	ex e1, e2;
	
	e1 = (1+x)*pow((2+x),2)*pow((3+x),3)*pow((4+x),4);
	e2 = sqrfree(expand(e1),lst(x));
	if (e1 != e2) {
		clog << "sqrfree(expand(" << e1 << ")) erroneously returned "
		     << e2 << endl;
		++result;
	}
	
	e1 = (x+y)*pow((x+2*y),2)*pow((x+3*y),3)*pow((x+4*y),4);
	e2 = sqrfree(expand(e1));
	if (e1 != e2) {
		clog << "sqrfree(expand(" << e1 << ")) erroneously returned "
		     << e2 << endl;
		++result;
	}
	e2 = sqrfree(expand(e1),lst(x));
	if (e1 != e2) {
		clog << "sqrfree(expand(" << e1 << "),[x]) erroneously returned "
		     << e2 << endl;
		++result;
	}
	e2 = sqrfree(expand(e1),lst(y));
	if (e1 != e2) {
		clog << "sqrfree(expand(" << e1 << "),[y]) erroneously returned "
		     << e2 << endl;
		++result;
	}
	e2 = sqrfree(expand(e1),lst(x,y));
	if (e1 != e2) {
		clog << "sqrfree(expand(" << e1 << "),[x,y]) erroneously returned "
		     << e2 << endl;
		++result;
	}
	
	return result;
}

/* Arithmetic Operators should behave just as one expects from built-in types.
 * When somebody screws up the operators this routine will most probably fail
 * to compile.  Unfortunately we can only test the stuff that is allowed, not
 * what is forbidden (e.g. e1+e2 = 42) since that must not compile.  :-(   */
static unsigned exam_operator_semantics()
{
	unsigned result = 0;
	ex e1, e2;
	int i1, i2;
	
	// Assignment should not return const ex though it may be obfuscated:
	e1 = 7; e2 = 4;
	i1 = 7; i2 = 4;
	(e1 = e2) = 2;
	(i1 = i2) = 2;
	if (e1!=i1 || e2!=i2) {
		clog << "Semantics of ex::operator=() screwed." << endl;
		++result;
	}
	(e1 += e2) = 2;
	(i1 += i2) = 2;
	if (e1!=i1 || e2!=i2) {
		clog << "Semantics of ex::operator=() screwed." << endl;
		++result;
	}
	(e1 -= e2) = 2;
	(i1 -= i2) = 2;
	if (e1!=i1 || e2!=i2) {
		clog << "Semantics of ex::operator=() screwed." << endl;
		++result;
	}
	
	// Prefix/postfix increment/decrement behaviour:
	e1 = 7; e2 = 4;
	i1 = 7; i2 = 4;
	e1 = (--e2 = 2)++;
	i1 = (--i2 = 2)++;
	if (e1!=i1 || e2!=i2) {
		clog << "Semantics of increment/decrement operators screwed." << endl;
		++result;
	}
	e1 = (++e2 = 2)--;
	i1 = (++i2 = 2)--;
	if (e1!=i1 || e2!=i2) {
		clog << "Semantics of increment/decrement operators screwed." << endl;
		++result;
	}
	
	// prefix increment/decrement must return an lvalue (contrary to postfix):
	e1 = 7; e2 = 4;
	i1 = 7; i2 = 4;
	--++----e1;  ++(++++++++(++++e2));
	--++----i1;  ++(++++++++(++++i2));
	if (e1!=i1 || e2!=i2) {
		clog << "Semantics of prefix increment/decrement operators screwed." << endl;
		++result;
	}
	
	// This one has a good chance of detecting problems in self-assignment:
	// (which incidentally was severely broken from version 0.7.3 to 0.8.2).
	ex selfprobe = numeric("65536");
	selfprobe = selfprobe;
	if (!is_exactly_a<numeric>(selfprobe)) {
		clog << "ex (of numeric) after self-assignment became " << selfprobe << endl;
		++result;
	}
	
	return result;
}

/* This checks whether subs() works as intended in some special cases. */
static unsigned exam_subs()
{
	unsigned result = 0;
	symbol x("x");
	ex e1, e2;

	// This used to fail in GiNaC 1.0.5 because it first substituted
	// x+1 -> (x-1)+1 -> x, and then substituted again x -> x-1, giving
	// the wrong result
	e1 = x+1;
	e2 = e1.subs(x == x-1);
	if (!e2.is_equal(x)) {
		clog << "(x+1).subs(x==x-1) erroneously returned " << e2 << " instead of x" << endl;
		++result;
	}

	// And this used to fail in GiNaC 1.5.8 because it first substituted
	// exp(x) -> exp(log(x)) -> x, and then substitued again x -> log(x)
	e1 = exp(x);
	e2 = e1.subs(x == log(x));
	if (!e2.is_equal(x)) {
		clog << "exp(x).subs(x==log(x)) erroneously returned " << e2 << " instead of x" << endl;
		++result;
	}

	e1 = sin(1+sin(x));
	e2 = e1.subs(sin(wild()) == cos(wild()));
	if (!e2.is_equal(cos(1+cos(x)))) {
		clog << "sin(1+sin(x)).subs(sin($1)==cos($1)) erroneously returned " << e2 << " instead of cos(1+cos(x))" << endl;
		++result;
	}

	return result;
}

/* Joris van der Hoeven (he of TeXmacs fame) is a funny guy.  He has his own
 * ideas what a symbolic system should do.  Let's make sure we won't disappoint
 * him some day.  Incidentally, this seems to always have worked. */
static unsigned exam_joris()
{
	unsigned result = 0;
	symbol x("x");

	ex e = expand(pow(x, x-1) * x);
	if (e != pow(x, x)) {
		clog << "x^(x-1)*x did not expand to x^x.  Please call Joris!" << endl;
		++result;
	}

	return result;
}

/* Test Chris Dams' algebraic substitutions. */
static unsigned exam_subs_algebraic()
{
	unsigned result = 0;
	symbol x("x"), y("y");

	ex e = ex(x*x*x*y*y).subs(x*y==2, subs_options::algebraic);
	if (e != 4*x) {
		clog << "(x^3*y^2).subs(x*y==2,subs_options::algebraic) erroneously returned " << e << endl;
		++result;
	}

	e = ex(x*x*x*x*x).subs(x*x==y, subs_options::algebraic);
	if (e != y*y*x) {
		clog << "x^5.subs(x^2==y,subs_options::algebraic) erroneously returned " << e << endl;
		++result;
	}

	e=x*x*y;
	if (!e.has(x*y, has_options::algebraic))
	{	clog << "(x^2*y).has(x*y, has_options::algebraic) erroneously returned false." << endl;
		++result;
	}

	if (e.has(x*y*y, has_options::algebraic))
	{	clog << "(x^2*y).has(x*y*y, has_options::algebraic) erroneously returned true." << endl;
		++result;
	}

	e=x*x*x*y;
	if (!e.has(x*x, has_options::algebraic))
	{	clog << "(x^3*y).has(x*x, has_options::algebraic) erroneously returned false." << endl;
		++result;
	}

	if (e.has(y*y, has_options::algebraic))
	{	clog << "(x^3*y).has(y*y, has_options::algebraic) erroneously returned true." << endl;
		++result;
	}

	return result;
}

unsigned exam_misc()
{
	unsigned result = 0;
	
	cout << "examining miscellaneous other things" << flush;
	
	result += exam_expand_subs();  cout << '.' << flush;
	result += exam_expand_subs2();  cout << '.' << flush;
	result += exam_expand_power(); cout << '.' << flush;
	result += exam_sqrfree(); cout << '.' << flush;
	result += exam_operator_semantics(); cout << '.' << flush;
	result += exam_subs(); cout << '.' << flush;
	result += exam_joris(); cout << '.' << flush;
	result += exam_subs_algebraic(); cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_misc();
}
