/** @file exam_numeric.cpp
 *
 *  These exams creates some numbers and check the result of several boolean
 *  tests on these numbers like is_integer() etc... */

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
#include <sstream>
using namespace std;

/* Simple and maybe somewhat pointless consistency tests of assorted tests and
 * conversions. */
static unsigned exam_numeric1()
{
	unsigned result = 0;
	numeric test_int1(42);
	numeric test_int2(5);
	numeric test_rat1 = test_int1; test_rat1 /= test_int2;
	test_rat1 = -test_rat1;		 // -42/5
	numeric test_crat = test_rat1+I*test_int2;  // 5*I-42/5
	symbol a("a");
	ex e1, e2;
	
	if (!test_int1.is_integer()) {
		clog << test_int1
		     << " erroneously not recognized as integer" << endl;
		++result;
	}
	if (!test_int1.is_rational()) {
		clog << test_int1
		     << " erroneously not recognized as rational" << endl;
		++result;
	}
	
	if (!test_rat1.is_rational()) {
		clog << test_rat1
		     << " erroneously not recognized as rational" << endl;
		++result;
	}
	if (test_rat1.is_integer()) {
		clog << test_rat1
			 << " erroneously recognized as integer" << endl;
		++result;
	}
	
	if (!test_crat.is_crational()) {
		clog << test_crat
		     << " erroneously not recognized as complex rational" << endl;
		++result;
	}
	if (test_crat.info(info_flags::nonnegative)) {
		clog << test_crat
		     << " erroneously recognized as non-negative number" << endl;
		++result;
	}
	
	int i = numeric(1984).to_int();
	if (i-1984) {
		clog << "conversion of " << i
		     << " from numeric to int failed" << endl;
		++result;
	}
	
	e1 = test_int1;
	if (!e1.info(info_flags::posint)) {
		clog << "expression " << e1
		     << " erroneously not recognized as positive integer" << endl;
		++result;
	}
	
	e2 = test_int1 + a;
	if (e2.info(info_flags::integer)) {
		clog << "expression " << e2
		     << " erroneously recognized as integer" << endl;
		++result;
	}
	
	// The next two were two actual bugs in CLN till June, 12, 1999:
	test_rat1 = numeric(3)/numeric(2);
	test_rat1 += test_rat1;
	if (!test_rat1.is_integer()) {
		clog << "3/2 + 3/2 erroneously not integer 3 but instead "
		     << test_rat1 << endl;
		++result;
	}
	test_rat1 = numeric(3)/numeric(2);
	numeric test_rat2 = test_rat1 + numeric(1);  // 5/2
	test_rat2 -= test_rat1;  // 1
	if (!test_rat2.is_integer()) {
		clog << "5/2 - 3/2 erroneously not integer 1 but instead "
		     << test_rat2 << endl;
		++result;
	}
	
	return result;
}

/* We had some fun with a bug in CLN that caused it to loop forever when
 * calculating expt(a,b) if b is a rational and a a nonnegative integer.
 * Implementing a workaround sadly introduced another bug on May 28th 1999
 * that was fixed on May 31st.  The workaround turned out to be stupid and
 * the original bug in CLN was finally killed on September 2nd. */
static unsigned exam_numeric2()
{
	unsigned result = 0;
	
	ex zero = numeric(0);
	ex two = numeric(2);
	ex three = numeric(3);
	
	// The hang in this code was the reason for the original workaround
	if (pow(two,two/three)==42) {
		clog << "pow(2,2/3) erroneously returned 42" << endl;
		++result;  // cannot happen
	}
	
	// Actually, this used to raise a FPE after introducing the workaround
	if (two*zero!=zero) {
		clog << "2*0 erroneously returned " << two*zero << endl;
		++result;
	}
	
	// And this returned a cl_F due to the implicit call of numeric::power()
	ex six = two*three;
	if (!six.info(info_flags::integer)) {
		clog << "2*3 erroneously returned the non-integer " << six << endl;
		++result;
	}
	
	// The fix in the workaround left a whole which was fixed hours later...
	ex another_zero = pow(zero,numeric(1)/numeric(2));
	if (!another_zero.is_zero()) {
		clog << "pow(0,1/2) erroneously returned" << another_zero << endl;
		++result;
	}
	
	return result;
}

/* Assorted tests to ensure some crucial functions behave exactly as specified
 * in the documentation. */
static unsigned exam_numeric3()
{
	unsigned result = 0;
	numeric calc_rem, calc_quo;
	numeric a, b;
	
	// check if irem(a, b) and irem(a, b, q) really behave like Maple's 
	// irem(a, b) and irem(a, b, 'q') as advertised in our documentation.
	// These overloaded routines indeed need to be checked separately since
	// internally they might be doing something completely different:
	a = 23; b = 4; calc_rem = irem(a, b);
	if (calc_rem != 3) {
		clog << "irem(" << a << "," << b << ") erroneously returned "
		     << calc_rem << endl;
		++result;
	}
	a = 23; b = -4; calc_rem = irem(a, b);
	if (calc_rem != 3) {
		clog << "irem(" << a << "," << b << ") erroneously returned "
			 << calc_rem << endl;
		++result;
	}
	a = -23; b = 4; calc_rem = irem(a, b);
	if (calc_rem != -3) {
		clog << "irem(" << a << "," << b << ") erroneously returned "
		     << calc_rem << endl;
		++result;
	}
	a = -23; b = -4; calc_rem = irem(a, b);
	if (calc_rem != -3) {
		clog << "irem(" << a << "," << b << ") erroneously returned "
		     << calc_rem << endl;
		++result;
	}
	// and now the overloaded irem(a,b,q):
	a = 23; b = 4; calc_rem = irem(a, b, calc_quo);
	if (calc_rem != 3 || calc_quo != 5) {
		clog << "irem(" << a << "," << b << ",q) erroneously returned "
		     << calc_rem << " with q=" << calc_quo << endl;
		++result;
	}
	a = 23; b = -4; calc_rem = irem(a, b, calc_quo);
	if (calc_rem != 3 || calc_quo != -5) {
		clog << "irem(" << a << "," << b << ",q) erroneously returned "
		     << calc_rem << " with q=" << calc_quo << endl;
		++result;
	}
	a = -23; b = 4; calc_rem = irem(a, b, calc_quo);
	if (calc_rem != -3 || calc_quo != -5) {
		clog << "irem(" << a << "," << b << ",q) erroneously returned "
		     << calc_rem << " with q=" << calc_quo << endl;
		++result;
	}
	a = -23; b = -4; calc_rem = irem(a, b, calc_quo);
	if (calc_rem != -3 || calc_quo != 5) {
		clog << "irem(" << a << "," << b << ",q) erroneously returned "
		     << calc_rem << " with q=" << calc_quo << endl;
		++result;
	}
	// check if iquo(a, b) and iquo(a, b, r) really behave like Maple's 
	// iquo(a, b) and iquo(a, b, 'r') as advertised in our documentation.
	// These overloaded routines indeed need to be checked separately since
	// internally they might be doing something completely different:
	a = 23; b = 4; calc_quo = iquo(a, b);
	if (calc_quo != 5) {
		clog << "iquo(" << a << "," << b << ") erroneously returned "
		     << calc_quo << endl;
		++result;
	}
	a = 23; b = -4; calc_quo = iquo(a, b);
	if (calc_quo != -5) {
		clog << "iquo(" << a << "," << b << ") erroneously returned "
		     << calc_quo << endl;
		++result;
	}
	a = -23; b = 4; calc_quo = iquo(a, b);
	if (calc_quo != -5) {
		clog << "iquo(" << a << "," << b << ") erroneously returned "
		     << calc_quo << endl;
		++result;
	}
	a = -23; b = -4; calc_quo = iquo(a, b);
	if (calc_quo != 5) {
		clog << "iquo(" << a << "," << b << ") erroneously returned "
		     << calc_quo << endl;
		++result;
	}
	// and now the overloaded iquo(a,b,r):
	a = 23; b = 4; calc_quo = iquo(a, b, calc_rem);
	if (calc_quo != 5 || calc_rem != 3) {
		clog << "iquo(" << a << "," << b << ",r) erroneously returned "
		     << calc_quo << " with r=" << calc_rem << endl;
		++result;
	}
	a = 23; b = -4; calc_quo = iquo(a, b, calc_rem);
	if (calc_quo != -5 || calc_rem != 3) {
		clog << "iquo(" << a << "," << b << ",r) erroneously returned "
		     << calc_quo << " with r=" << calc_rem << endl;
		++result;
	}
	a = -23; b = 4; calc_quo = iquo(a, b, calc_rem);
	if (calc_quo != -5 || calc_rem != -3) {
		clog << "iquo(" << a << "," << b << ",r) erroneously returned "
		     << calc_quo << " with r=" << calc_rem << endl;
		++result;
	}
	a = -23; b = -4; calc_quo = iquo(a, b, calc_rem);
	if (calc_quo != 5 || calc_rem != -3) {
		clog << "iquo(" << a << "," << b << ",r) erroneously returned "
		     << calc_quo << " with r=" << calc_rem << endl;
		++result;
	}
	
	return result;
}

/* Now we perform some less trivial checks about several functions which should
 * return exact numbers if possible. */
static unsigned exam_numeric4()
{
	unsigned result = 0;
	bool passed;
	
	// square roots of squares of integers:
	passed = true;
	for (int i=0; i<42; ++i)
		if (!sqrt(numeric(i*i)).is_integer())
			passed = false;
	if (!passed) {
		clog << "One or more square roots of squares of integers did not return exact integers" << endl;
		++result;
	}
	
	// square roots of squares of rationals:
	passed = true;
	for (int num=0; num<41; ++num)
		for (int den=1; den<42; ++den)
			if (!sqrt(numeric(num*num)/numeric(den*den)).is_rational())
				passed = false;
	if (!passed) {
		clog << "One or more square roots of squares of rationals did not return exact integers" << endl;
		++result;
	}
	
	return result;
}

/* This test examines that simplifications of the form 5^(3/2) -> 5*5^(1/2)
 * are carried out properly. */
static unsigned exam_numeric5()
{
	unsigned result = 0;
	
	// A variation of one of Ramanujan's wonderful identities must be
	// verifiable with very primitive means:
	ex e1 = pow(1 + pow(3,numeric(1,5)) - pow(3,numeric(2,5)),3);
	ex e2 = expand(e1 - 10 + 5*pow(3,numeric(3,5)));
	if (!e2.is_zero()) {
		clog << "expand((1+3^(1/5)-3^(2/5))^3-10+5*3^(3/5)) returned "
		     << e2 << " instead of 0." << endl;
		++result;
	}
	
	return result;
}

/* This test checks whether the numeric output/parsing routines are
   consistent. */
static unsigned exam_numeric6()
{
	unsigned result = 0;

	symbol sym("sym");
	vector<ex> test_numbers;
	test_numbers.push_back(numeric(0));			// zero
	test_numbers.push_back(numeric(1));			// one
	test_numbers.push_back(numeric(-1));		// minus one
	test_numbers.push_back(numeric(42));		// positive integer
	test_numbers.push_back(numeric(-42));		// negative integer
	test_numbers.push_back(numeric(14,3));		// positive rational
	test_numbers.push_back(numeric(-14,3));		// negative rational
	test_numbers.push_back(numeric(3.141));		// positive decimal
	test_numbers.push_back(numeric(-3.141));	// negative decimal
	test_numbers.push_back(numeric(0.1974));	// positive decimal, leading zero
	test_numbers.push_back(numeric(-0.1974));	// negative decimal, leading zero
	test_numbers.push_back(sym);				// symbol

	for (vector<ex>::const_iterator br=test_numbers.begin(); br<test_numbers.end(); ++br) {
		for (vector<ex>::const_iterator bi=test_numbers.begin(); bi<test_numbers.end(); ++bi) {

			for (vector<ex>::const_iterator er=test_numbers.begin(); er<test_numbers.end(); ++er) {
				for (vector<ex>::const_iterator ei=test_numbers.begin(); ei<test_numbers.end(); ++ei) {

					// Construct expression, don't test invalid ones
					ex base = (*br) + (*bi)*I, exponent = (*er) + (*ei)*I, x;
					try {
						x = pow(base, exponent);
					} catch (...) {
						continue;
					}

					// Print to string
					std::ostringstream s;
					s << x;

					// Read back expression from string
					string x_as_string = s.str();
					ex x_again(x_as_string, lst(sym));

					// They should be equal
					if (!x_again.is_equal(x)) {
						clog << x << " was read back as " << x_again << endl;
						++result;
					}
				}
			}
		}
	}

	return result;
}

unsigned exam_numeric()
{
	unsigned result = 0;
	
	cout << "examining consistency of numeric types" << flush;
	
	result += exam_numeric1();  cout << '.' << flush;
	result += exam_numeric2();  cout << '.' << flush;
	result += exam_numeric3();  cout << '.' << flush;
	result += exam_numeric4();  cout << '.' << flush;
	result += exam_numeric5();  cout << '.' << flush;
	result += exam_numeric6();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_numeric();
}
