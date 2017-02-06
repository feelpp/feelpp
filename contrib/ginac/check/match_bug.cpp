/** @file match_bug.cpp
 *
 *  Check for bug in GiNaC::ex::match() described here:
 *  http://www.ginac.de/pipermail/ginac-devel/2006-April/000942.html */

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
#include "error_report.h"

#include <iostream>
using namespace GiNaC;

/*
 * basic::match(lst&) used to have an obscure side effect: repl_lst
 * could be modified even if the match failed! Although this "feature"
 * was documented it happened to be very confusing *even for GiNaC
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

/*
 * expairseq::match() should not have any side effects if the match failed.
 */
static void expairseq_failed_match_no_side_effect(int count)
{
	for (int i = 0; i < count; ++i) {
		exmap repls;
		symbol t("t"), A("A");
		ex e = pow(t, 2)*exp(t*A);
		ex pattern = pow(t, wild(0))*exp(wild(1))*A;
		bool matched = e.match(pattern, repls);
		cbug_on(matched, "unexpected match: " << e << " vs " << pattern);
		cbug_on(repls.size(), "failed match has side effects");
	}
}

/*
 * exp(a)*sin(x) + exp(b)*sin(y) used to fail to match
 * exp(a)*sin($0) + exp(b)*sin($1). The failure was not deterministic.
 *
 * The first attempted submatch is sin(y)*exp(b) with sin($0)*exp(a).
 * It fails but $0 == y gets assigned due to a spurious side effect.
 * Next submatch is sin(x)*exp(a) with sin($0)*exp(a) (the same pattern
 * as in the first submatch). This one fails because of (incorrect)
 * $0 == y assignment.
 *
 * Note: due to the unstable term ordering the sequence of submatches
 * might be different and the match might succeed (as it should), hence
 * we repeat the test several times.
 */
static void expairseq_match_false_negative(int count)
{
	for (int i = 0; i < count; ++i) {
		symbol a("a"), b("b"), x("x"), y("y");
		ex e = exp(a)*sin(x) + exp(b)*sin(y);
		ex pattern = exp(a)*sin(wild(0)) + exp(b)*sin(wild(1));
		cbug_on(!e.match(pattern), "match failed: " << e << "did not"
			"match " << pattern);
	}
}

int main(int argc, char** argv)
{
	const int repetitions = 100;
	std::cout << "checking for historical bugs in match()... " << std::flush;
	failed_match_have_side_effects();
	match_false_negative();
	expairseq_failed_match_no_side_effect(repetitions);
	expairseq_match_false_negative(repetitions);
	std::cout << "not found. ";
	return 0;
}
