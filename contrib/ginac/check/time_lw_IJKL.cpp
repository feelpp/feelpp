/** @file time_lw_IJKL.cpp
 *
 *  Tests I, J, K and L from the paper "Comparison of Polynomial-Oriented CAS"
 *  by Robert H. Lewis and Michael Wester. */

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

static unsigned test(unsigned n)
{
	unsigned result = 0;
	timer cartier;
	char name = (n==40?'I':(n==70?'K':'?'));
	
	cout << "timing Lewis-Wester test " << name
	     << " (invert rank " << n << " Hilbert)" << flush;
	
	// Create a rank n Hilbert matrix:
	matrix H(n,n);
	for (unsigned r=0; r<n; ++r)
		for (unsigned c=0; c<n; ++c)
			H.set(r,c,numeric(1,r+c+1));
	// invert it:
	cartier.start();
	matrix Hinv(n,n);
	Hinv = H.inverse();
	cout << ". passed ";
	cout << cartier.read() << 's' << endl;
	
	// check result:
	name = (n==40?'J':(n==70?'L':'?'));
	
	cout << "timing Lewis-Wester test " << name
	     << " (check rank " << n << " Hilbert)" << flush;
	
	cartier.reset();
	matrix identity = H.mul(Hinv);
	bool correct = true;
	for (unsigned r=0; r<n; ++r)
		for (unsigned c=0; c<n; ++c) {
			if (r==c) {
				if (identity(r,c)!=1)
					correct = false;
			} else {
				if (identity(r,c)!=0)
					correct = false;
			}
		}
	if (!correct)
		++result;
	cout << cartier.read() << 's' << endl;
	return result;
}

unsigned time_lw_IJKL()
{
	unsigned result = 0;
	
	// Tests I and J:
	result += test(40);
	// Tests K and L:
	result += test(70);
	
	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_lw_IJKL();
}
