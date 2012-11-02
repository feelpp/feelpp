/** @file numeric_archive.cpp
 *
 *  Check for a bug in numeric::archive
 *
 *  numeric::archive used to fail if the real part of a complex number
 *  is a rational number and the imaginary part is a floating point one. */

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

#include <algorithm>
#include <cln/cln.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>
using namespace cln;

struct archive_unarchive_check
{
	cl_N operator()(const cl_N& n) const
	{
		ex e = numeric(n);
		archive ar;
		ar.archive_ex(e, "test");
		lst l;
		ex check = ar.unarchive_ex(l, "test");
		if (!check.is_equal(e)) {
			std::ostringstream s;
			s << __FILE__ << ':' << __LINE__ << ": expected: " << e << ", got " << check;
			throw std::logic_error(s.str());
		}
		return n;
	}
};

int main(int argc, char** argv)
{
	const cl_I one(1);
	std::cout << "checking if numeric::archive handles complex numbers properly" << std::endl;
	const cl_R three_fp = cl_float(3.0, default_float_format);
	std::vector<cl_N> numbers;
	numbers.push_back(complex(one, three_fp));
	numbers.push_back(complex(three_fp, one));
	numbers.push_back(complex(three_fp, three_fp));
	numbers.push_back(complex(one, one));
	std::for_each(numbers.begin(), numbers.end(), archive_unarchive_check());
	return 0;
}
