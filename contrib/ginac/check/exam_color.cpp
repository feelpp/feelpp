/** @file exam_color.cpp
 *
 *  Here we test GiNaC's color objects (su(3) Lie algebra). */

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

static unsigned check_equal(const ex &e1, const ex &e2)
{
	ex e = e1 - e2;
	if (!e.is_zero()) {
		clog << e1 << "-" << e2 << " erroneously returned "
		     << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned check_equal_simplify(const ex &e1, const ex &e2)
{
	ex e = simplify_indexed(e1) - e2;
	if (!e.is_zero()) {
		clog << "simplify_indexed(" << e1 << ")-" << e2 << " erroneously returned "
		     << e << " instead of 0" << endl;
		return 1;
	}
	return 0;
}

static unsigned color_check1()
{
	// checks general identities and contractions of the structure constants

	unsigned result = 0;

	idx a(symbol("a"), 8), b(symbol("b"), 8), c(symbol("c"), 8), d(symbol("d"), 8);

	result += check_equal(color_d(a, c, a), 0);
	result += check_equal_simplify(color_d(a, b, c) * color_d(b, d, c), numeric(5,3) * delta_tensor(a, d));
	result += check_equal_simplify(color_d(idx(5, 8), b, c) * color_d(b, idx(5, 8), c), numeric(5,3));
	result += check_equal_simplify(color_d(a, b, c) * color_d(b, c, a), numeric(40,3));
	result += check_equal_simplify(color_d(a, b, c) * color_f(b, d, c), 0);
	result += check_equal_simplify(color_d(a, b, c) * color_f(b, c, a), 0);
	result += check_equal_simplify(color_f(a, b, c) * color_f(b, c, a), 24);
	result += check_equal_simplify(color_f(a, b, c) * color_f(b, d, c), -3 * delta_tensor(a, d));
	result += check_equal_simplify(color_h(a, b, c) * color_h(a, b, c), numeric(-32,3));
	result += check_equal_simplify(color_h(a, b, c) * color_h(b, a, c), numeric(112,3));

    ex e = color_h(a, b, c) * color_h(a, b, c);
    ex sum = 0;
	for (int i=1; i<9; i++)
	    for (int j=1; j<9; j++)
	        for (int k=1; k<9; k++)
                sum += e.subs(lst(a == i, b == j, c == k));
	if (!sum.is_equal(numeric(-32,3))) {
		clog << "numeric contraction of " << e << " erroneously returned "
		     << sum << " instead of -32/3" << endl;
		result++;
	}

	return result;
}

static unsigned color_check2()
{
	// checks general identities and contractions of the generators

	unsigned result = 0;

	idx a(symbol("a"), 8), b(symbol("b"), 8), c(symbol("c"), 8), k(symbol("k"), 8);
	ex e;

	e = color_T(k) * color_T(k);
	result += check_equal_simplify(e, 4 * color_ONE() / 3);
	e = color_T(k) * color_T(a) * color_T(k);
	result += check_equal_simplify(e, -color_T(a) / 6);
	e = color_T(k) * color_T(a) * color_T(b) *  color_T(k);
	result += check_equal_simplify(e, delta_tensor(a, b) * color_ONE() / 4 - color_T(a) * color_T(b) / 6);
	e = color_T(k) * color_T(a) * color_T(b) *  color_T(c) * color_T(k);
	result += check_equal_simplify(e, (color_h(a, b, c) * color_ONE() / 8).expand() - color_T(a) * color_T(b) * color_T(c) / 6);
	e = color_T(a) * color_T(b) * color_T(a) *  color_T(b);
	result += check_equal_simplify(e, -2 * color_ONE() / 9);
	e = color_T(a) * color_T(b) * color_T(b) *  color_T(a);
	result += check_equal_simplify(e, 16 * color_ONE() / 9);
	e = color_T(a) * color_T(b) * color_T(c) * color_T(c) * color_T(b) *  color_T(a);
	result += check_equal_simplify(e, 64 * color_ONE() / 27);
	e = color_T(a) * color_T(b) * color_T(c) * color_T(k) * color_T(a) * color_T(k) *  color_T(c) * color_T(b);
	result += check_equal_simplify(e, -color_ONE() / 162);

	return result;
}

static unsigned color_check3()
{
	// checks traces

	unsigned result = 0;

	idx a(symbol("a"), 8), b(symbol("b"), 8), c(symbol("c"), 8);
	ex e;

	e = color_ONE();
	result += check_equal(color_trace(e), 3);
	e = color_T(a);
	result += check_equal(color_trace(e), 0);
	e = color_T(a) * color_T(b);
	result += check_equal(color_trace(e), delta_tensor(a, b) / 2);
	e = color_T(a) * color_T(b) * color_T(c);
	result += check_equal(color_trace(e), color_h(a, b, c) / 4);

	e = color_ONE(0) * color_ONE(1) / 9;
	result += check_equal(color_trace(e, 0), color_ONE(1) / 3);
	result += check_equal(color_trace(e, 1), color_ONE(0) / 3);
	result += check_equal(color_trace(e, 2), e);
	result += check_equal(color_trace(e, lst(0, 1)), 1);

	e = color_T(a, 0) * color_T(a, 1) * color_T(b, 0) * color_T(b, 1);
	result += check_equal_simplify(color_trace(e, 0), 2 * color_ONE(1) / 3);
	result += check_equal_simplify(color_trace(e, 1), 2 * color_ONE(0) / 3);
	result += check_equal_simplify(color_trace(e, 2), e);
	result += check_equal_simplify(color_trace(e, lst(0, 1)), 2);

	return result;
}

unsigned exam_color()
{
	unsigned result = 0;
	
	cout << "examining color objects" << flush;

	result += color_check1();  cout << '.' << flush;
	result += color_check2();  cout << '.' << flush;
	result += color_check3();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_color();
}
