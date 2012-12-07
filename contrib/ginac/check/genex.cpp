/** @file genex.cpp
 *
 *  Provides some routines for generating expressions that are later used as 
 *  input in the consistency checks. */

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

#include <cstdlib>
using namespace std;

/* Create a dense univariate random polynomial in x.
 * (of the form 9 - 22*a - 17*a^2 + 14*a^3 + 7*a^4 + 7a^5 if degree==5) */
const ex
dense_univariate_poly(const symbol & x, unsigned degree)
{
	ex unipoly;
	
	for (unsigned i=0; i<=degree; ++i)
		unipoly += numeric((rand()-RAND_MAX/2))*pow(x,i);
	
	return unipoly;
}

/* Create a dense bivariate random polynomial in x1 and x2.
 * (of the form 9 + 52*x1 - 27*x1^2 + 84*x2 + 7*x2^2 - 12*x1*x2 if degree==2)
 */
const ex
dense_bivariate_poly(const symbol & x1, const symbol & x2, unsigned degree)
{
	ex bipoly;
	
	for (unsigned i1=0; i1<=degree; ++i1)
		for (unsigned i2=0; i2<=degree-i1; ++i2)
			bipoly += numeric((rand()-RAND_MAX/2))*pow(x1,i1)*pow(x2,i2);
	
	return bipoly;
}

/* Chose a randum symbol or number from the argument list. */
const ex
random_symbol(const symbol & x,
			  const symbol & y,
			  const symbol & z,
			  bool rational = true,
			  bool complex = false)
{
	ex e;
	switch (abs(rand()) % 4) {
		case 0:
			e = x;
			break;
		case 1:
			e = y;
			break;
		case 2:
			e = z;
			break;
		case 3: {
			int c1;
			do { c1 = rand()%20 - 10; } while (!c1);
			int c2;
			do { c2 = rand()%20 - 10; } while (!c2);
			if (!rational)
				c2 = 1;
			e = numeric(c1, c2);
			if (complex && !(rand()%5))
				e = e*I;
			break;
		}
	}
	return e;
}

/* Create a sparse random tree in three symbols. */
const ex
sparse_tree(const symbol & x,
			const symbol & y,
			const symbol & z,
			int level,
			bool trig = false,	// true includes trigonomatric functions
			bool rational = true, // false excludes coefficients in Q
			bool complex = false) // true includes complex numbers
{
	if (level == 0)
		return random_symbol(x,y,z,rational,complex);
	switch (abs(rand()) % 10) {
		case 0:
		case 1:
		case 2:
		case 3:
			return add(sparse_tree(x,y,z,level-1, trig, rational),
					   sparse_tree(x,y,z,level-1, trig, rational));
		case 4:
		case 5:
		case 6:
			return mul(sparse_tree(x,y,z,level-1, trig, rational),
					   sparse_tree(x,y,z,level-1, trig, rational));
		case 7:
		case 8: {
			ex powbase;
			do {
				powbase = sparse_tree(x,y,z,level-1, trig, rational);
			} while (powbase.is_zero());
			return pow(powbase, abs(rand() % 4));
			break;
		}
		case 9:
			if (trig) {
				switch (abs(rand()) % 4) {
					case 0:
						return sin(sparse_tree(x,y,z,level-1, trig, rational));
					case 1:
						return cos(sparse_tree(x,y,z,level-1, trig, rational));
					case 2:
						return exp(sparse_tree(x,y,z,level-1, trig, rational));
					case 3: {
						ex logex;
						do {
							ex logarg;
							do {
								logarg = sparse_tree(x,y,z,level-1, trig, rational);
							} while (logarg.is_zero());
							// Keep the evaluator from accidentally plugging an
							// unwanted I in the tree:
							if (!complex && logarg.info(info_flags::negative))
								logarg = -logarg;
							logex = log(logarg);
						} while (logex.is_zero());
						return logex;
						break;
					}
				}
			} else
				return random_symbol(x,y,z,rational,complex);
	}

	return 0;
}
