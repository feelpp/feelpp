/** @file time_lw_N.cpp
 *
 *  Test N from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
 *  Lewis and Michael Wester (also known as the smaller version of the first
 *  Fermat-test). */

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
#include <vector>
using namespace std;

static const bool do_test = false;  // set to true in order to run this beast

static unsigned test()
{
	symbol p11("p11"), p12("p12"), p21("p21"), p22("p22");
	symbol a12("a12"), a21("a21"), a22("a22");
	symbol n11("n11"), n22("n22");
	symbol g("g");
	symbol q1("q1"), q2("q2"), q3("q3"), q4("q4");
	
	ex ss1 = ex("(4*g*a22^3-g*a12*a21*a22^2-4*n22*a21*a22^2+4*n11*a21*a22^2+7*g*a12*a22^2+4*g^2*a22^2-4*n11*n22*a22^2+4*n11*a22^2+n22*a12*a21^2*a22-n11*a12*a21^2*a22-4*g*a21^2*a22-g*a12^2*a21*a22-5*g^2*a12*a21*a22+5*n11*n22*a12*a21*a22-7*n22*a12*a21*a22+2*n11*a12*a21*a22-4*g*a21*a22+3*g*a12^2*a22+3*g^2*a12*a22-3*n11*n22*a12*a22+3*n11*a12*a22+g*a12*a21^3+g^2*a12^2*a21^2-n11*n22*a12^2*a21^2+n22*a12^2*a21^2-2*g*a12*a21^2-3*g^2*a12^2*a21+3*n11*n22*a12^2*a21-3*n22*a12^2*a21-3*g*a12*a21)/(3*g*a12*a21*a22^2-3*n22*a21*a22^2+g*a12*a22^2-n22*a22^2+3*n11*a12*a21^2*a22-3*g*a21^2*a22+5*g*a12^2*a21*a22-5*n22*a12*a21*a22+4*n11*a12*a21*a22-4*g*a21*a22+g*a12^2*a22-n22*a12*a22+n11*a12*a22-g*a22+2*n11*a12^2*a21^2-2*g*a12*a21^2+2*g*a12^3*a21-2*n22*a12^2*a21+2*n11*a12^2*a21-2*g*a12*a21)",lst(g,a12,a21,a22,n11,n22));
	
	ex ss2 = ex("(4*g*a12*a22^2-4*n22*a22^2+4*a22^2-g*a12^2*a21*a22+n22*a12*a21*a22+4*n11*a12*a21*a22-5*a12*a21*a22-4*g*a21*a22+3*g*a12^2*a22-3*n22*a12*a22+3*a12*a22-n11*a12^2*a21^2+a12^2*a21^2+g*a12*a21^2+3*n11*a12^2*a21-3*a12^2*a21-3*g*a12*a21)/(2*g*a12*a22^2-2*n22*a22^2+g*a12^2*a21*a22-n22*a12*a21*a22+2*n11*a12*a21*a22-2*g*a21*a22+2*g*a12^2*a22-2*n22*a12*a22+2*n11*a12*a22-2*g*a22+n11*a12^2*a21^2-g*a12*a21^2+g*a12^3*a21-n22*a12^2*a21+n11*a12^2*a21-g*a12*a21)",lst(g,a12,a21,a22,n11,n22));
	
	ex ss3 = ex("(4*p21*a22^3-p21*a12*a21*a22^2-4*p22*a21*a22^2+4*p11*a21*a22^2+7*p21*a12*a22^2-4*p11*p22*a22^2+4*p12*p21*a22^2+4*p11*a22^2+p22*a12*a21^2*a22-p11*a12*a21^2*a22-4*p12*a21^2*a22-p21*a12^2*a21*a22+5*p11*p22*a12*a21*a22-7*p22*a12*a21*a22-5*p12*p21*a12*a21*a22+2*p11*a12*a21*a22-4*p12*a21*a22+3*p21*a12^2*a22-3*p11*p22*a12*a22+3*p12*p21*a12*a22+3*p11*a12*a22+p12*a12*a21^3-p11*p22*a12^2*a21^2+p22*a12^2*a21^2+p12*p21*a12^2*a21^2-2*p12*a12*a21^2+3*p11*p22*a12^2*a21-3*p22*a12^2*a21-3*p12*p21*a12^2*a21-3*p12*a12*a21)/(3*p21*a12*a21*a22^2-3*p22*a21*a22^2+p21*a12*a22^2-p22*a22^2+3*p11*a12*a21^2*a22-3*p12*a21^2*a22+5*p21*a12^2*a21*a22-5*p22*a12*a21*a22+4*p11*a12*a21*a22-4*p12*a21*a22+p21*a12^2*a22-p22*a12*a22+p11*a12*a22-p12*a22+2*p11*a12^2*a21^2-2*p12*a12*a21^2+2*p21*a12^3*a21-2*p22*a12^2*a21+2*p11*a12^2*a21-2*p12*a12*a21)",lst(a12,a21,a22,p11,p12,p21,p22));
	
	ex ss4 = ex("(4*p21*a12*a22^2-4*p22*a22^2+4*a22^2-p21*a12^2*a21*a22+p22*a12*a21*a22+4*p11*a12*a21*a22-5*a12*a21*a22-4*p12*a21*a22+3*p21*a12^2*a22-3*p22*a12*a22+3*a12*a22-p11*a12^2*a21^2+a12^2*a21^2+p12*a12*a21^2+3*p11*a12^2*a21-3*a12^2*a21-3*p12*a12*a21)/(2*p21*a12*a22^2-2*p22*a22^2+p21*a12^2*a21*a22-p22*a12*a21*a22+2*p11*a12*a21*a22-2*p12*a21*a22+2*p21*a12^2*a22-2*p22*a12*a22+2*p11*a12*a22-2*p12*a22+p11*a12^2*a21^2-p12*a12*a21^2+p21*a12^3*a21-p22*a12^2*a21+p11*a12^2*a21-p12*a12*a21)",lst(p11,p12,p21,p22,a12,a21,a22));
	
	ex res1 = ex("p11*p22*q1^2*q4^2-p12*p21*q1^2*q4^2-n22*p11*p22*q1*q4^2+2*n11*p11*p22*q1*q4^2-p11*p22*q1*q4^2+n22*p12*p21*q1*q4^2-2*n11*p12*p21*q1*q4^2+p12*p21*q1*q4^2+2*g^2*p11*p22*q4^2-2*n11*n22*p11*p22*q4^2+2*n22*p11*p22*q4^2+2*n11*p11*p22*q4^2-2*p11*p22*q4^2-2*g^2*p12*p21*q4^2+2*n11*n22*p12*p21*q4^2-2*n22*p12*p21*q4^2-2*n11*p12*p21*q4^2+2*p12*p21*q4^2-n11*p22*q1*q2*q3*q4+g*p21*q1*q2*q3*q4+g*p12*q1*q2*q3*q4-n22*p11*q1*q2*q3*q4-g^2*p22*q2*q3*q4+n11*n22*p22*q2*q3*q4-n11*p22*q2*q3*q4-2*g*p21*q2*q3*q4+g*p12*q2*q3*q4+2*g^2*p11*q2*q3*q4-2*n11*n22*p11*q2*q3*q4+2*n22*p11*q2*q3*q4-n11*p22*q1*q3*q4+p22*q1*q3*q4+g*p21*q1*q3*q4-2*g*p12*q1*q3*q4+2*n22*p11*q1*q3*q4-2*p11*q1*q3*q4-g^2*p22*q3*q4+n11*n22*p22*q3*q4-n22*p22*q3*q4-n11*p22*q3*q4+p22*q3*q4-4*g^2*p11*q3*q4+4*n11*n22*p11*q3*q4-4*n22*p11*q3*q4-4*n11*p11*q3*q4+4*p11*q3*q4+n22*p11*p22*q1*q2*q4-2*n11*p11*p22*q1*q2*q4+2*n11*p22*q1*q2*q4-n22*p12*p21*q1*q2*q4+2*n11*p12*p21*q1*q2*q4+g*p21*q1*q2*q4-2*g*p12*q1*q2*q4-n22*p11*q1*q2*q4-4*g^2*p11*p22*q2*q4+4*n11*n22*p11*p22*q2*q4-2*n22*p11*p22*q2*q4-2*n11*p11*p22*q2*q4+2*g^2*p22*q2*q4-2*n11*n22*p22*q2*q4+2*n11*p22*q2*q4+4*g^2*p12*p21*q2*q4-4*n11*n22*p12*p21*q2*q4+2*n22*p12*p21*q2*q4+2*n11*p12*p21*q2*q4-2*g*p21*q2*q4-2*g*p12*q2*q4+2*g^2*p11*q2*q4-2*n11*n22*p11*q2*q4+2*n22*p11*q2*q4-p11*p22*q1^2*q4-p22*q1^2*q4+p12*p21*q1^2*q4+2*p11*q1^2*q4-n22*p11*p22*q1*q4-4*n11*p11*p22*q1*q4+5*p11*p22*q1*q4+n22*p22*q1*q4-p22*q1*q4+n22*p12*p21*q1*q4+4*n11*p12*p21*q1*q4-5*p12*p21*q1*q4+g*p21*q1*q4+4*g*p12*q1*q4+4*n11*p11*q1*q4-4*p11*q1*q4-g^2*q2^2*q3^2+n11*n22*q2^2*q3^2+g^2*q2*q3^2-n11*n22*q2*q3^2-n22*q2*q3^2+2*n11*q2*q3^2+2*g^2*q3^2-2*n11*n22*q3^2+2*n22*q3^2+2*n11*q3^2-2*q3^2+g^2*p22*q2^2*q3-n11*n22*p22*q2^2*q3-2*g^2*p11*q2^2*q3+2*n11*n22*p11*q2^2*q3+g^2*q2^2*q3-n11*n22*q2^2*q3+2*n11*p22*q1*q2*q3-2*g*p21*q1*q2*q3+g*p12*q1*q2*q3-n22*p11*q1*q2*q3+n22*q1*q2*q3-2*n11*q1*q2*q3+g^2*p22*q2*q3-n11*n22*p22*q2*q3+n22*p22*q2*q3+4*g*p21*q2*q3+g*p12*q2*q3+4*g^2*p11*q2*q3-4*n11*n22*p11*q2*q3+4*n11*p11*q2*q3-5*g^2*q2*q3+5*n11*n22*q2*q3-n22*q2*q3-4*n11*q2*q3+2*n11*p22*q1*q3-2*p22*q1*q3-2*g*p21*q1*q3-2*g*p12*q1*q3+2*n22*p11*q1*q3-2*p11*q1*q3-2*n22*q1*q3-2*n11*q1*q3+4*q1*q3+2*g^2*p11*p22*q2^2-2*n11*n22*p11*p22*q2^2-2*g^2*p22*q2^2+2*n11*n22*p22*q2^2-2*g^2*p12*p21*q2^2+2*n11*n22*p12*p21*q2^2-2*g^2*p11*q2^2+2*n11*n22*p11*q2^2+2*g^2*q2^2-2*n11*n22*q2^2+n22*p11*p22*q1*q2+4*n11*p11*p22*q1*q2-n22*p22*q1*q2-4*n11*p22*q1*q2-n22*p12*p21*q1*q2-4*n11*p12*p21*q1*q2-n22*p11*q1*q2-4*n11*p11*q1*q2+n22*q1*q2+4*n11*q1*q2-2*p11*p22*q1^2+2*p22*q1^2+2*p12*p21*q1^2+2*p11*q1^2-2*q1^2",lst(p11,p12,p21,p22,n11,n22,g,q1,q2,q3,q4));
	ex result = res1.subs(lst(q1==ss1, q2==ss2, q3==ss3, q4==ss4));
	ex normalresult = normal(result);
	if (!normalresult.is_zero()) {
		clog << "Normalization should have returned 0." << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_N()
{
	unsigned result = 0;
	unsigned count = 0;
	timer tag_heuer;
	double time = .0;
	
	cout << "timing Lewis-Wester test N (poly at rational fcns)" << flush;
	
	if (do_test) {
		tag_heuer.start();
		// correct for very small times:
		do {
			result = test();
			++count;
		} while ((time=tag_heuer.read())<0.1 && !result);
		cout << '.' << flush;
		cout << time/count << 's' << endl;
	} else {
		cout << " disabled" << endl;
	}
	
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_N();
}
