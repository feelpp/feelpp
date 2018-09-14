/** @file exam_misc.cpp
 *
 */

/*
 *  GiNaC Copyright (C) 1999-2016 Johannes Gutenberg University Mainz, Germany
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

/* Exam real/imaginary part of polynomials. */
static unsigned exam_polynomials()
{
	realsymbol a("a"), b("b");
	ex e = pow(a + I*b,3).expand() + I;

	if (e.real_part() != pow(a,3)-3*a*pow(b,2) ||
	    e.imag_part() != 1+3*pow(a,2)*b-pow(b,3)) {
		clog << "real / imaginary part miscomputed" << endl;
		return 1;
	}
	return 0;
}

/* Exam symbolic expansion of nested expression. */
static unsigned exam_monster()
{
	// This little monster is inspired by Sage's Symbench R1.
	// It is much more aggressive that the original and covers more code.
	struct {  // C++ doesn't have nested functions...
		ex operator()(const ex & z) {
			return sqrt(ex(1)/3) * pow(z, 11) - I / pow(2*z, 3);
		}
	} f;
	ex monster = f(f(f(f(I/2))));  // grows exponentially with number of nestings..
	ex r = real_part(monster);
	ex i = imag_part(monster);

	// Check with precomputed result
	double r_eps = ex_to<numeric>(evalf(r)).to_double() - 0.2000570104163233;
	double i_eps = ex_to<numeric>(evalf(i)).to_double() - 0.5284320312415462;
	if (abs(r_eps) > 1e-9  ||  abs(i_eps) > 1e-9) {
		clog << "iterated function was miscomputed" << endl;
		return 1;
	}
	return 0;
}

unsigned exam_real_imag()
{
	unsigned result = 0;

	cout << "examining real/imaginary part separation" << flush;

	result += exam_polynomials(); cout << '.' << flush;
	result += exam_monster(); cout << '.' << flush;

	return result;
}

int main(int argc, char** argv)
{
	return exam_real_imag();
}
