/** @file exam_matrices.cpp
 *
 *  Here we examine manipulations on GiNaC's symbolic matrices. */

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
#include <stdexcept>
using namespace std;

static unsigned matrix_determinants()
{
	unsigned result = 0;
	ex det;
	matrix m1(1,1), m2(2,2), m3(3,3), m4(4,4);
	symbol a("a"), b("b"), c("c");
	symbol d("d"), e("e"), f("f");
	symbol g("g"), h("h"), i("i");
	
	// check symbolic trivial matrix determinant
	m1 = matrix{{a}};
	det = m1.determinant();
	if (det != a) {
		clog << "determinant of 1x1 matrix " << m1
		     << " erroneously returned " << det << endl;
		++result;
	}
	
	// check generic dense symbolic 2x2 matrix determinant
	m2 = matrix{{a, b},
	            {c, d}};
	det = m2.determinant();
	if (det != (a*d-b*c)) {
		clog << "determinant of 2x2 matrix " << m2
		     << " erroneously returned " << det << endl;
		++result;
	}
	
	// check generic dense symbolic 3x3 matrix determinant
	m3 = matrix{{a, b, c},
	            {d, e, f},
	            {g, h, i}};
	det = m3.determinant();
	if (det != (a*e*i - a*f*h - d*b*i + d*c*h + g*b*f - g*c*e)) {
		clog << "determinant of 3x3 matrix " << m3
		     << " erroneously returned " << det << endl;
		++result;
	}
	
	// check dense numeric 3x3 matrix determinant
	m3 = matrix{{0, -1,  3},
	            {3, -2,  2},
	            {3,  4, -2}};
	det = m3.determinant();
	if (det != 42) {
		clog << "determinant of 3x3 matrix " << m3
		     << " erroneously returned " << det << endl;
		++result;
	}
	
	// check dense symbolic 2x2 matrix determinant
	m2 = matrix{{a/(a-b), 1},
	            {b/(a-b), 1}};
	det = m2.determinant();
	if (det != 1) {
		if (det.normal() == 1)  // only half wrong
			clog << "determinant of 2x2 matrix " << m2
			     << " was returned unnormalized as " << det << endl;
		else  // totally wrong
			clog << "determinant of 2x2 matrix " << m2
			     << " erroneously returned " << det << endl;
		++result;
	}
	
	// check sparse symbolic 4x4 matrix determinant
	m4.set(0,1,a).set(1,0,b).set(3,2,c).set(2,3,d);
	det = m4.determinant();
	if (det != a*b*c*d) {
		clog << "determinant of 4x4 matrix " << m4
		     << " erroneously returned " << det << endl;
		++result;
	}
	
	// check characteristic polynomial
	m3 = matrix{{a, -2,   2},
		    {3, a-1,  2},
		    {3,  4,  a-3}};
	ex p = m3.charpoly(a);
	if (p != 0) {
		clog << "charpoly of 3x3 matrix " << m3
		     << " erroneously returned " << p << endl;
		++result;
	}
	
	return result;
}

static unsigned matrix_invert1()
{
	unsigned result = 0;
	matrix m(1,1);
	symbol a("a");
	
	m.set(0,0,a);
	matrix m_i = m.inverse();
	
	if (m_i(0,0) != pow(a,-1)) {
		clog << "inversion of 1x1 matrix " << m
		     << " erroneously returned " << m_i << endl;
		++result;
	}
	
	return result;
}

static unsigned matrix_invert2()
{
	unsigned result = 0;
	symbol a("a"), b("b"), c("c"), d("d");
	matrix m = {{a, b},
	            {c, d}};
	matrix m_i = m.inverse();
	ex det = m.determinant();
	
	if ((normal(m_i(0,0)*det) != d) ||
	    (normal(m_i(0,1)*det) != -b) ||
	    (normal(m_i(1,0)*det) != -c) ||
	    (normal(m_i(1,1)*det) != a)) {
		clog << "inversion of 2x2 matrix " << m
		     << " erroneously returned " << m_i << endl;
		++result;
	}
	
	return result;
}

static unsigned matrix_invert3()
{
	unsigned result = 0;
	symbol a("a"), b("b"), c("c");
	symbol d("d"), e("e"), f("f");
	symbol g("g"), h("h"), i("i");
	matrix m = {{a, b, c},
	            {d, e, f},
	            {g, h, i}};
	matrix m_i = m.inverse();
	ex det = m.determinant();
	
	if ((normal(m_i(0,0)*det) != (e*i-f*h)) ||
	    (normal(m_i(0,1)*det) != (c*h-b*i)) ||
	    (normal(m_i(0,2)*det) != (b*f-c*e)) ||
	    (normal(m_i(1,0)*det) != (f*g-d*i)) ||
	    (normal(m_i(1,1)*det) != (a*i-c*g)) ||
	    (normal(m_i(1,2)*det) != (c*d-a*f)) ||
	    (normal(m_i(2,0)*det) != (d*h-e*g)) ||
	    (normal(m_i(2,1)*det) != (b*g-a*h)) ||
	    (normal(m_i(2,2)*det) != (a*e-b*d))) {
		clog << "inversion of 3x3 matrix " << m
		     << " erroneously returned " << m_i << endl;
		++result;
	}
	
	return result;
}

static unsigned matrix_solve2()
{
	// check the solution of the multiple system A*X = B:
	//       [ 1  2 -1 ] [ x0 y0 ]   [ 4 0 ]
	//       [ 1  4 -2 ]*[ x1 y1 ] = [ 7 0 ]
	//       [ a -2  2 ] [ x2 y2 ]   [ a 4 ]
	unsigned result = 0;
	symbol a("a");
	symbol x0("x0"), x1("x1"), x2("x2");
	symbol y0("y0"), y1("y1"), y2("y2");
	matrix A = {{1,  2, -1},
	            {1,  4, -2},
	            {a, -2,  2}};
	matrix B = {{4, 0},
	            {7, 0},
	            {a, 4}};
	matrix X = {{x0 ,y0},
	            {x1, y1},
	            {x2, y2}};
	matrix cmp = {{1, 0},
	              {3, 2},
	              {3, 4}};
	matrix sol(A.solve(X, B));
	if (cmp != sol) {
		clog << "Solving " << A << " * " << X << " == " << B << endl
		     << "erroneously returned " << sol << endl;
		result = 1;
	}

	return result;
}

static unsigned matrix_evalm()
{
	unsigned result = 0;

	matrix S {{1, 2},
	          {3, 4}};
	matrix T {{1, 1},
	          {2, -1}};
	matrix R {{27, 14},
	          {36, 26}};

	ex e = ((S + T) * (S + 2*T));
	ex f = e.evalm();
	if (!f.is_equal(R)) {
		clog << "Evaluating " << e << " erroneously returned " << f << " instead of " << R << endl;
		result++;
	}

	return result;
}

static unsigned matrix_rank()
{
	unsigned result = 0;
	symbol x("x"), y("y");
	matrix m(3,3);

	// the zero matrix always has rank 0
	if (m.rank() != 0) {
		clog << "The rank of " << m << " was not computed correctly." << endl;
		++result;
	}

	// a trivial rank one example
	m = {{1, 0, 0},
	     {2, 0, 0},
	     {3, 0, 0}};
	if (m.rank() != 1) {
		clog << "The rank of " << m << " was not computed correctly." << endl;
		++result;
	}

	// an example from Maple's help with rank two
	m = {{x,   1,  0},
	     {0,   0,  1},
	     {x*y, y,  1}};
	if (m.rank() != 2) {
		clog << "The rank of " << m << " was not computed correctly." << endl;
		++result;
	}

	// the 3x3 unit matrix has rank 3
	m = ex_to<matrix>(unit_matrix(3,3));
	if (m.rank() != 3) {
		clog << "The rank of " << m << " was not computed correctly." << endl;
		++result;
	}

	return result;	
}

static unsigned matrix_misc()
{
	unsigned result = 0;
	symbol a("a"), b("b"), c("c"), d("d"), e("e"), f("f");
	matrix m1 = {{a, b},
		     {c, d}};
	ex tr = trace(m1);
	
	// check a simple trace
	if (tr.compare(a+d)) {
		clog << "trace of 2x2 matrix " << m1
		     << " erroneously returned " << tr << endl;
		++result;
	}
	
	// and two simple transpositions
	matrix m2 = transpose(m1);
	if (m2(0,0) != a || m2(0,1) != c || m2(1,0) != b || m2(1,1) != d) {
		clog << "transpose of 2x2 matrix " << m1
			 << " erroneously returned " << m2 << endl;
		++result;
	}
	matrix m3 = {{a, b},
		     {c, d},
		     {e, f}};
	if (transpose(transpose(m3)) != m3) {
		clog << "transposing 3x2 matrix " << m3 << " twice"
		     << " erroneously returned " << transpose(transpose(m3)) << endl;
		++result;
	}
	
	// produce a runtime-error by inverting a singular matrix and catch it
	matrix m4(2,2);
	matrix m5;
	bool caught = false;
	try {
		m5 = inverse(m4);
	} catch (std::runtime_error err) {
		caught = true;
	}
	if (!caught) {
		cerr << "singular 2x2 matrix " << m4
		     << " erroneously inverted to " << m5 << endl;
		++result;
	}
	
	return result;
}

unsigned exam_matrices()
{
	unsigned result = 0;
	
	cout << "examining symbolic matrix manipulations" << flush;
	
	result += matrix_determinants();  cout << '.' << flush;
	result += matrix_invert1();  cout << '.' << flush;
	result += matrix_invert2();  cout << '.' << flush;
	result += matrix_invert3();  cout << '.' << flush;
	result += matrix_solve2();  cout << '.' << flush;
	result += matrix_evalm();  cout << "." << flush;
	result += matrix_rank();  cout << "." << flush;
	result += matrix_misc();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_matrices();
}
