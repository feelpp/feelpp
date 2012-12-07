/** @file match_bug.cpp
 *
 *  Check for bug in GiNaC::ex::match() described here:
 *  http://www.ginac.de/pipermail/ginac-devel/2006-April/000942.html */

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
#include "error_report.h"

#include <iostream>
using namespace GiNaC;

/*
 * basic::match(lst&) used to have an obscure side effect: repl_lst
 * could be modified even if the match failed! Although this "feature"
 * was documented it happend to be very confusing *even for GiNaC
 * developers*, see 
 * http://www.ginac.de/pipermail/ginac-devel/2006-April/000942.html
 *
 * It was fixed in 192ed7390b7b2b705ad100e3db0a92eedd2b20ad. Let's make
 * sure it will be never re-added:
 */
static void failed_match_have_side_effects()
{
	symbol x("x");
	ex e = pow(x, 5);
	ex pattern = pow(wild(0), -1);
	// obviously e does NOT match the pattern
	exmap repls;
	bool match_p = e.match(pattern, repls);
	cbug_on(match_p, "match(" << e << ", " << pattern << ") says \"Yes\"");
	cbug_on(repls.size() != 0,
		"failed match have side effects: repls = " << repls);
}

/*
 * As a consequence of the bug described above pattern matching can wrongly
 * fail. In particular, x^5*y^(-1) fails to match ($0)^(-1)*x^($2).
 *
 * The first thing that is attempted to match is x^5 with $0^(-1). This match
 * will fail. However repl_lst will contain $0 == x as a side effect. This
 * repl_lst will prevent the match of y^(-1) to ($0)^(-1) to succeed.
 *
 * This issue was worked around by 73f0ce4cf8d91f073f35a45443f5fbe886921c5c.
 * Now we have a real fix (192ed7390b7b2b705ad100e3db0a92eedd2b20ad), but
 * let's add a check.
 */
static void match_false_negative()
{
	symbol x("x"), y("y");
	ex e = pow(x, 5)*pow(y, -1);
	ex pattern = pow(wild(0), -1)*pow(x, wild(2));
	exmap repls;
	bool match_p = e.match(pattern, repls);
	cbug_on(!match_p, "false negative: " << e << " did not match "
			<< pattern);
}

int main(int argc, char** argv)
{
	std::cout << "checking for historical bugs in match()... " << std::flush;
	failed_match_have_side_effects();
	match_false_negative();
	std::cout << "not found. ";
	return 0;
}
