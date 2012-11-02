/** @file time_lw_M1.cpp
 *
 *  Test M1 from the paper "Comparison of Polynomial-Oriented CAS" by Robert H.
 *  Lewis and Michael Wester. */

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
using namespace std;

static unsigned test()
{
	// Determinant of a sparse matrix that comes up in graph theory:
	symbol x1("x1"), x2("x2"), x3("x3"), x4("x4"), x5("x5");
	ex w[26][11] = {
		{ 1,  1,  1,  7, x4, 12, x3, 17, x2, 22, x1},
		{ 2,  2,  1,  8, x4, 13, x3, 18, x2, 23, x1},
		{ 3,  3,  1,  9, x4, 14, x3, 19, x2, 24, x1},
		{ 4,  4,  1, 10, x4, 15, x3, 20, x2, 25, x1},
		{ 5,  5,  1, 26,  1,  1,  0,  1,  0,  1, 0 },
		{ 6,  2, x5,  6,  1, 12, x3, 17, x2, 22, x1},
		{ 7,  3, x5,  7,  1, 13, x3, 18, x2, 23, x1},
		{ 8,  4, x5,  8,  1, 14, x3, 19, x2, 24, x1},
		{ 9,  5, x5,  9,  1, 15, x3, 20, x2, 25, x1},
		{10, 10,  1, 26,  1,  1,  0,  1,  0,  1, 0 },
		{11,  2, x5,  7, x4, 11,  1, 17, x2, 22, x1},
		{12,  3, x5,  8, x4, 12,  1, 18, x2, 23, x1},
		{13,  4, x5,  9, x4, 13,  1, 19, x2, 24, x1},
		{14,  5, x5, 10, x4, 14,  1, 20, x2, 25, x1},
		{15, 15,  1, 26,  1,  1,  0,  1,  0,  1, 0 },
		{16,  2, x5,  7, x4, 12, x3, 16,  1, 22, x1},
		{17,  3, x5,  8, x4, 13, x3, 17,  1, 23, x1},
		{18,  4, x5,  9, x4, 14, x3, 18,  1, 24, x1},
		{19,  5, x5, 10, x4, 15, x3, 19,  1, 25, x1},
		{20, 20,  1, 26,  1,  1,  0,  1,  0,  1, 0 },
		{21,  2, x5,  7, x4, 12, x3, 17, x2, 21, 1 },
		{22,  3, x5,  8, x4, 13, x3, 18, x2, 22, 1 },
		{23,  4, x5,  9, x4, 14, x3, 19, x2, 23, 1 },
		{24,  5, x5, 10, x4, 15, x3, 20, x2, 24, 1 },
		{25, 25,  1, 26,  1,  1,  0,  1,  0,  1, 0 },
		{26,  1, x5,  6, x4, 11, x3, 16, x2, 21, x1}
	};
	matrix m(26,26);
	for (unsigned r=0; r<26; ++r) {
		for (unsigned c=0; c<5; ++c) {
			m.set(r,
			      unsigned(ex_to<numeric>(w[r][2*c+1]).to_int()-1),
			      w[r][2*c+2]);
		}
	}
	ex det = m.determinant();
	// The result should have been:
	ex cmp("-12*x2^2*x5^2*x4-12*x1*x5^2*x3^2-x5^3*x4^2-12*x1*x5^2*x4^2-12*x2*x5^2*x4^2-12*x3*x5^2*x4^2-x4^3*x5^2-36*x3*x1*x5^2*x4-36*x3*x1*x4^2*x5-36*x3*x2*x5^2*x4-36*x3*x2*x4^2*x5-2*x5^3*x4*x2-12*x3^2*x5^2*x4-12*x3^2*x4^2*x5-2*x5^3*x4*x3-2*x4^3*x5*x3-12*x1*x5^2*x2^2-36*x1*x5*x3^2*x4-36*x2*x5*x3^2*x4-x3^3*x5^2-x3^3*x4^2-2*x3^3*x5*x4-12*x2^2*x4^2*x5-12*x2*x5^2*x3^2-12*x2*x4^2*x3^2-12*x1*x4^2*x3^2-x3^2*x5^3-x3^2*x4^3-2*x4^3*x5*x2-2*x3*x5^3*x2-2*x3*x4^3*x2-2*x3^3*x5*x2-2*x3^3*x4*x2-2*x2^3*x5*x4-2*x2^3*x5*x3-2*x2^3*x4*x3-36*x2^2*x5*x4*x3-36*x2*x1*x5^2*x4-36*x2*x1*x4^2*x5-120*x2*x1*x5*x4*x3-36*x2*x1*x5^2*x3-36*x2*x1*x4^2*x3-36*x2*x1*x3^2*x5-36*x2*x1*x3^2*x4-12*x2^2*x5^2*x3-12*x2^2*x4^2*x3-12*x2^2*x3^2*x5-12*x2^2*x3^2*x4-2*x1^3*x4*x3-2*x1^3*x4*x2-2*x1^3*x3*x2-2*x1^3*x5*x2-36*x1^2*x5*x4*x3-36*x2*x1^2*x5*x4-36*x2*x3*x1^2*x5-36*x2*x3*x1^2*x4-x1^3*x5^2-x1^3*x4^2-x1^3*x3^2-x1^3*x2^2-x2^2*x5^3-x2^2*x4^3-x2^2*x3^3-12*x1*x4^2*x2^2-12*x1*x3^2*x2^2-12*x1^2*x5^2*x4-12*x1^2*x4^2*x5-12*x1^2*x5^2*x3-12*x1^2*x4^2*x3-12*x1^2*x3^2*x5-12*x1^2*x3^2*x4-12*x1^2*x5^2*x2-12*x1^2*x4^2*x2-12*x1^2*x3^2*x2-12*x1^2*x2^2*x5-12*x1^2*x2^2*x4-12*x1^2*x2^2*x3-2*x5^3*x4*x1-2*x4^3*x5*x1-2*x3*x5^3*x1-2*x3*x4^3*x1-2*x3^3*x5*x1-2*x3^3*x4*x1-2*x2*x5^3*x1-2*x2*x4^3*x1-2*x2*x3^3*x1-2*x2^3*x5*x1-2*x2^3*x4*x1-2*x2^3*x3*x1-2*x1^3*x5*x4-2*x1^3*x5*x3-36*x1*x5*x2^2*x4-36*x1*x5*x2^2*x3-36*x1*x4*x2^2*x3-x1^2*x5^3-x1^2*x4^3-x1^2*x3^3-x2^3*x5^2-x2^3*x4^2-x2^3*x3^2-x1^2*x2^3",lst(x1,x2,x3,x4,x5));
	if (det!=cmp) {
		clog << "The determinant was miscalculated" << endl;
		return 1;
	}
	return 0;
}

unsigned time_lw_M1()
{
	unsigned result = 0;
	unsigned count = 0;
	timer rolex;
	double time = .0;
	
	cout << "timing Lewis-Wester test M1 (26x26 sparse, det)" << flush;
	
	rolex.start();
	// correct for very small times:
	do {
		result = test();
		++count;
	} while ((time=rolex.read())<0.1 && !result);
	cout << '.' << flush;
	cout << time/count << 's' << endl;
	
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_M1();
}
