/** @file parser_memleak.cpp
 *
 *  This small program exhibits the memory leak in the ginac_yylex().
 *  Run it as
 *
 *  valgrind --leak-check=yes  ./parser_memleak
 *
 *  or simply
 *
 *  ulimit -v `expr 64 \* 1024` ./parser_memleak
 */

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

#include <ginac/ginac.h>
using namespace GiNaC;

#include <iostream>
#include <stdexcept>
using namespace std;

int main(int argc, char** argv) {
	const symbol x("x"), y("y");
	const lst syms(x, y);
	// parser-generated symbol => memory leak.
	static const char* str[] = { "x^2+2*x*y + cos(x)", "Li2(x/y) + log(y/x)" };
	
	// depends on the amount of the available VM, compiler options, etc.
	const unsigned N_max = 500000;
	unsigned N=0;
	ex e;
	try {
		for (; N < N_max; N++) {
			e = ex(str[N & 1], syms);
		}
	} catch (std::bad_alloc) {
		cerr << "N = " << N << endl;
		return 1;
	}
	return 0;
}
