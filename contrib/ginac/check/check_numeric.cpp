/** @file check_numeric.cpp
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

#include <cstdlib> // for rand()
#include <iostream>
using namespace std;

/* Simple and maybe somewhat pointless consistency tests of assorted tests and
 * conversions. */
static unsigned check_numeric1()
{
	unsigned result = 0;
	bool errorflag = false;
	int re_q, im_q;
	
	// Check some numerator and denominator calculations:
	for (int rep=0; rep<200; ++rep) {
		do { re_q = rand(); } while (re_q == 0);
		do { im_q = rand(); } while (im_q == 0);
		numeric r(rand()-RAND_MAX/2, re_q);
		numeric i(rand()-RAND_MAX/2, im_q);
		numeric z = r + I*i;
		numeric p = numer(z);
		numeric q = denom(z);
		numeric res = p/q;
		if (res != z) {
			clog << z << " erroneously transformed into " 
			     << p << "/" << q << " by numer() and denom()" << endl;
			errorflag = true;
		}
	}
	if (errorflag)
		++result;
	
	return result;
}

static unsigned check_numeric2()
{
	unsigned result = 0;
	bool errorflag = false;
	int i_num, i_den;
	
	// Check non-nested radicals (n/d)^(m/n) in ex wrapper class:
	for (int i=0; i<200; ++i) {
		for (int j=2; j<13; ++j) {
			// construct an exponent 1/j...
			numeric nm(1,j);
			nm += numeric(int(20.0*rand()/(RAND_MAX+1.0))-10);
			// ...a numerator...
			do {
				i_num = rand();
			} while (i_num<=0);
			numeric num(i_num);
			// ...and a denominator.
			do {
				i_den = (rand())/100;
			} while (i_den<=0);
			numeric den(i_den);
			// construct the radicals:
			ex radical = pow(ex(num)/ex(den),ex(nm));
			numeric floating = pow(num/den,nm);
			// test the result:
			if (is_a<numeric>(radical)) {
				// This is very improbable with decent random numbers but it
				// still can happen, so we better check if it is correct:
				if (pow(radical,inverse(nm))==num/den) {
					// Aha! We drew some lucky numbers. Nothing to see here...
				} else {
					clog << "(" << num << "/" << den << ")^(" << nm
						 << ") should have been a product, instead it's "
						 << radical << endl;
					errorflag = true;
				}
			}
			numeric ratio = abs(ex_to<numeric>(evalf(radical))/floating);
			if (ratio>1.0001 && ratio<0.9999) {
				clog << "(" << num << "/" << den << ")^(" << nm
				     << ") erroneously evaluated to " << radical;
				errorflag = true;
			}
		}
	}
	if (errorflag)
		++result;
	
	return result;
}

unsigned check_numeric()
{
	unsigned result = 0;
	
	cout << "checking consistency of numeric types" << flush;
	
	result += check_numeric1();  cout << '.' << flush;
	result += check_numeric2();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return check_numeric();
}
