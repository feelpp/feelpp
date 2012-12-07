/** @file upoly_io.cpp
 *
 *  Input/Output function for univariate polynomials. */

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

#include "upoly.h"
#include "upoly_io.h"

#include <cln/integer_io.h>
#include <iostream>
#include <string>

namespace GiNaC {

using std::ostream;
using std::string;

template<typename T> static void 
print(const T& p, ostream& os, const string& varname = string("x"))
{
	if (p.size() == 0)
		os << '0';

	bool seen_nonzero = false;

	for (std::size_t i = p.size(); i-- != 0;  ) {
		if (zerop(p[i])) {
			if (seen_nonzero)
				continue;
			os << "+ [WARNING: 0]*" << varname << "^" << i << "]";
			continue;
		}
		seen_nonzero = true;
		os << "+ (" << p[i] << ")";
		if (i != 0)
			os << "*" << varname;
		if (i > 1)
			os << '^' << i;
		os << " ";
	}
}

#define DEFINE_OPERATOR_OUT(type)                         \
std::ostream& operator<<(std::ostream& os, const type& p) \
{                                                         \
	print(p, os);                                         \
	return os;                                            \
}                                                         \
void dbgprint(const type& p)                              \
{                                                         \
	print(p, std::cerr);                                  \
}

DEFINE_OPERATOR_OUT(upoly);
DEFINE_OPERATOR_OUT(umodpoly);
#undef DEFINE_OPERATOR_OUT

} // namespace GiNaC
