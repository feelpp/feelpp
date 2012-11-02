/** @file exam_lsolve.cpp
 *
 *  These exams test solving small linear systems of symbolic equations. */

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

static unsigned exam_lsolve1()
{
	// A trivial example.
	unsigned result = 0;
	symbol x("x");
	ex eq, aux;
	
	eq = (3*x+5 == numeric(8));
	aux = lsolve(eq, x);
	if (aux != 1) {
		++result;
		clog << "solution of 3*x+5==8 erroneously returned "
		     << aux << endl;
	}
	
	return result;
}

static unsigned exam_lsolve2a()
{
	// An example from the Maple online help.
	unsigned result = 0;
	symbol a("a"), b("b"), x("x"), y("y");
	lst eqns, vars;
	ex sol;
	
	// Create the linear system [a*x+b*y==3,x-y==b]...
	eqns.append(a*x+b*y==3).append(x-y==b);
	// ...to be solved for [x,y]...
	vars.append(x).append(y);
	// ...and solve it:
	sol = lsolve(eqns, vars);
	ex sol_x = sol.op(0).rhs();  // rhs of solution for first variable (x)
	ex sol_y = sol.op(1).rhs();  // rhs of solution for second variable (y)
	
	// It should have returned [x==(3+b^2)/(a+b),y==(3-a*b)/(a+b)]
	if (!normal(sol_x - (3+pow(b,2))/(a+b)).is_zero() ||
		!normal(sol_y - (3-a*b)/(a+b)).is_zero()) {
		++result;
		clog << "solution of the system " << eqns << " for " << vars
		     << " erroneously returned " << sol << endl;
	}
	
	return result;
}

static unsigned exam_lsolve2b()
{
	// A boring example from Mathematica's online help.
	unsigned result = 0;
	symbol x("x"), y("y");
	lst eqns, vars;
	ex sol;
	
	// Create the linear system [3*x+y==7,2*x-5*y==8]...
	eqns.append(3*x+y==7).append(2*x-5*y==8);
	// ...to be solved for [x,y]...
	vars.append(x).append(y);
	// ...and solve it:
	sol = lsolve(eqns, vars);
	ex sol_x = sol.op(0).rhs();  // rhs of solution for first variable (x)
	ex sol_y = sol.op(1).rhs();  // rhs of solution for second variable (y)
	
	// It should have returned [x==43/17,y==-10/17]
	if ((sol_x != numeric(43,17)) ||
		(sol_y != numeric(-10,17))) {
		++result;
		clog << "solution of the system " << eqns << " for " << vars
		     << " erroneously returned " << sol << endl;
	}
	
	return result;
}

static unsigned exam_lsolve2c()
{
	// A more interesting example from the Maple online help.
	unsigned result = 0;
	symbol x("x"), y("y");
	lst eqns, vars;
	ex sol;
	
	// Create the linear system [I*x+y==1,I*x-y==2]...
	eqns.append(I*x+y==1).append(I*x-y==2);
	// ...to be solved for [x,y]...
	vars.append(x).append(y);
	// ...and solve it:
	sol = lsolve(eqns, vars);
	ex sol_x = sol.op(0).rhs();  // rhs of solution for first variable (x)
	ex sol_y = sol.op(1).rhs();  // rhs of solution for second variable (y)
	
	// It should have returned [x==-3/2*I,y==-1/2]
	if ((sol_x != numeric(-3,2)*I) ||
		(sol_y != numeric(-1,2))) {
		++result;
		clog << "solution of the system " << eqns << " for " << vars
		     << " erroneously returned " << sol << endl;
	}
	
	return result;
}

static unsigned exam_lsolve2S()
{
	// A degenerate example that went wrong in GiNaC 0.6.2.
	unsigned result = 0;
	symbol x("x"), y("y"), t("t");
	lst eqns, vars;
	ex sol;
	
	// Create the linear system [0*x+0*y==0,0*x+1*y==t]...
	eqns.append(0*x+0*y==0).append(0*x+1*y==t);
	// ...to be solved for [x,y]...
	vars.append(x).append(y);
	// ...and solve it:
	sol = lsolve(eqns, vars);
	ex sol_x = sol.op(0).rhs();  // rhs of solution for first variable (x)
	ex sol_y = sol.op(1).rhs();  // rhs of solution for second variable (y)
	
	// It should have returned [x==x,y==t]
	if ((sol_x != x) ||
		(sol_y != t)) {
		++result;
		clog << "solution of the system " << eqns << " for " << vars
		     << " erroneously returned " << sol << endl;
	}
	
	return result;
}

static unsigned exam_lsolve3S()
{
	// A degenerate example that went wrong while trying to improve elimination
	unsigned result = 0;
	symbol b("b"), c("c");
	symbol x("x"), y("y"), z("z");
	lst eqns, vars;
	ex sol;
	
	// Create the linear system [y+z==b,-y+z==c] with one additional row...
	eqns.append(ex(0)==ex(0)).append(b==z+y).append(c==z-y);
	// ...to be solved for [x,y,z]...
	vars.append(x).append(y).append(z);
	// ...and solve it:
	sol = lsolve(eqns, vars);
	ex sol_x = sol.op(0).rhs();  // rhs of solution for first variable (x)
	ex sol_y = sol.op(1).rhs();  // rhs of solution for second variable (y)
	ex sol_z = sol.op(2).rhs();  // rhs of solution for third variable (z)
	
	// It should have returned [x==x,y==t,]
	if ((sol_x != x) ||
		(sol_y != (b-c)/2) ||
		(sol_z != (b+c)/2)) {
		++result;
		clog << "solution of the system " << eqns << " for " << vars
		     << " erroneously returned " << sol << endl;
	}
	
	return result;
}

unsigned exam_lsolve()
{
	unsigned result = 0;
	
	cout << "examining linear solve" << flush;
	
	result += exam_lsolve1();  cout << '.' << flush;
	result += exam_lsolve2a();  cout << '.' << flush;
	result += exam_lsolve2b();  cout << '.' << flush;
	result += exam_lsolve2c();  cout << '.' << flush;
	result += exam_lsolve2S();  cout << '.' << flush;
	result += exam_lsolve3S();  cout << '.' << flush;
	
	return result;
}

int main(int argc, char** argv)
{
	return exam_lsolve();
}
