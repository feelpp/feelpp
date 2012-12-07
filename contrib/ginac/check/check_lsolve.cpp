/** @file check_lsolve.cpp
 *
 *  These test routines do some simple checks on solving linear systems of
 *  symbolic equations.  They are a well-tried resource for cross-checking
 *  the underlying symbolic manipulations. */

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
#include <sstream>
using namespace std;

extern const ex 
dense_univariate_poly(const symbol & x, unsigned degree);

static unsigned check_matrix_solve(unsigned m, unsigned n, unsigned p,
								   unsigned degree)
{
	const symbol a("a");
	matrix A(m,n);
	matrix B(m,p);
	// set the first min(m,n) rows of A and B
	for (unsigned ro=0; (ro<m)&&(ro<n); ++ro) {
		for (unsigned co=0; co<n; ++co)
			A.set(ro,co,dense_univariate_poly(a,degree));
		for (unsigned co=0; co<p; ++co)
			B.set(ro,co,dense_univariate_poly(a,degree));
	}
	// repeat excessive rows of A and B to avoid excessive construction of
	// overdetermined linear systems
	for (unsigned ro=n; ro<m; ++ro) {
		for (unsigned co=0; co<n; ++co)
			A.set(ro,co,A(ro-1,co));
		for (unsigned co=0; co<p; ++co)
			B.set(ro,co,B(ro-1,co));
	}
	// create a vector of n*p symbols all named "xrc" where r and c are ints
	vector<symbol> x;
	matrix X(n,p);
	for (unsigned i=0; i<n; ++i) {
		for (unsigned j=0; j<p; ++j) {
			ostringstream buf;
			buf << "x" << i << j << ends;
			x.push_back(symbol(buf.str()));
			X.set(i,j,x[p*i+j]);
		}
	}
	matrix sol(n,p);
	// Solve the system A*X==B:
	try {
		sol = A.solve(X, B);
	} catch (const exception & err) {  // catch runtime_error
		// Presumably, the coefficient matrix A was degenerate
		string errwhat = err.what();
		if (errwhat == "matrix::solve(): inconsistent linear system")
			return 0;
		else
			clog << "caught exception: " << errwhat << endl;
		throw;
	}
	
	// check the result with our original matrix:
	bool errorflag = false;
	for (unsigned ro=0; ro<m; ++ro) {
		for (unsigned pco=0; pco<p; ++pco) {
			ex e = 0;
			for (unsigned co=0; co<n; ++co)
			e += A(ro,co)*sol(co,pco);
			if (!(e-B(ro,pco)).normal().is_zero())
				errorflag = true;
		}
	}
	if (errorflag) {
		clog << "Our solve method claims that A*X==B, with matrices" << endl
		     << "A == " << A << endl
		     << "X == " << sol << endl
		     << "B == " << B << endl;
		return 1;
	}
	
	return 0;
}

static unsigned check_inifcns_lsolve(unsigned n)
{
	unsigned result = 0;
	
	for (int repetition=0; repetition<200; ++repetition) {
		// create two size n vectors of symbols, one for the coefficients
		// a[0],..,a[n], one for indeterminates x[0]..x[n]:
		vector<symbol> a;
		vector<symbol> x;
		for (unsigned i=0; i<n; ++i) {
			ostringstream buf;
			buf << i << ends;
			a.push_back(symbol(string("a")+buf.str()));
			x.push_back(symbol(string("x")+buf.str()));
		}
		lst eqns;  // equation list
		lst vars;  // variable list
		ex sol;	// solution
		// Create a random linear system...
		for (unsigned i=0; i<n; ++i) {
			ex lhs = rand()%201-100;
			ex rhs = rand()%201-100;
			for (unsigned j=0; j<n; ++j) {
				// ...with small coefficients to give degeneracy a chance...
				lhs += a[j]*(rand()%21-10);
				rhs += x[j]*(rand()%21-10);
			}
			eqns.append(lhs==rhs);
			vars.append(x[i]);
		}
		// ...solve it...
		sol = lsolve(eqns, vars);
		
		// ...and check the solution:
		if (sol.nops() == 0) {
			// no solution was found
			// is the coefficient matrix really, really, really degenerate?
			matrix coeffmat(n,n);
			for (unsigned ro=0; ro<n; ++ro)
				for (unsigned co=0; co<n; ++co)
					coeffmat.set(ro,co,eqns.op(co).rhs().coeff(a[co],1));
			if (!coeffmat.determinant().is_zero()) {
				++result;
				clog << "solution of the system " << eqns << " for " << vars
					 << " was not found" << endl;
			}
		} else {
			// insert the solution into rhs of out equations
			bool errorflag = false;
			for (unsigned i=0; i<n; ++i)
				if (eqns.op(i).rhs().subs(sol) != eqns.op(i).lhs())
					errorflag = true;
			if (errorflag) {
				++result;
				clog << "solution of the system " << eqns << " for " << vars
				     << " erroneously returned " << sol << endl;
			}
		}
	}
	
	return result;
}

unsigned check_lsolve()
{
	unsigned result = 0;
	
	cout << "checking linear solve" << flush;
	
	// solve some numeric linear systems
	for (unsigned n=1; n<14; ++n)
		result += check_matrix_solve(n, n, 1, 0);
	cout << '.' << flush;
	// solve some underdetermined numeric systems
	for (unsigned n=1; n<14; ++n)
		result += check_matrix_solve(n+1, n, 1, 0);
	cout << '.' << flush;
	// solve some overdetermined numeric systems
	for (unsigned n=1; n<14; ++n)
		result += check_matrix_solve(n, n+1, 1, 0);
	cout << '.' << flush;
	// solve some multiple numeric systems
	for (unsigned n=1; n<14; ++n)
		result += check_matrix_solve(n, n, n/3+1, 0);
	cout << '.' << flush;
	// solve some symbolic linear systems
	for (unsigned n=1; n<8; ++n)
		result += check_matrix_solve(n, n, 1, 2);
	cout << '.' << flush;
	
	// check lsolve, the wrapper function around matrix::solve()
	result += check_inifcns_lsolve(2);  cout << '.' << flush;
	result += check_inifcns_lsolve(3);  cout << '.' << flush;
	result += check_inifcns_lsolve(4);  cout << '.' << flush;
	result += check_inifcns_lsolve(5);  cout << '.' << flush;
	result += check_inifcns_lsolve(6);  cout << '.' << flush;
		
	return result;
}

int main(int argc, char** argv)
{
	return check_lsolve();
}
