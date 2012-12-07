/** @file test_runner.h
 *
 *  Utility functions for benchmarking. */

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

#ifndef GINAC_CHECK_TEST_RUNNER_H
#define GINAC_CHECK_TEST_RUNNER_H

#include "timer.h"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

template<typename T> static void 
run_benchmark(T& benchmark,
	      const unsigned ntests = 10,
	      const std::clock_t t_annoying = 15*CLOCKS_PER_SEC,
	      const std::clock_t t_min = CLOCKS_PER_SEC/100)
{
	const std::clock_t start(std::clock());
	std::clock_t elapsed(start);

	timer P36;
	unsigned n = 0;
	double t = 0;
	bool go_on = true;
	do {
		++n;
		P36.start();
		benchmark.run();
		t += P36.read();
		go_on = benchmark.check();
		if (!go_on)
			break;
		elapsed = std::clock() - start;
		if (elapsed > t_annoying)
			break;
	} while (n <= ntests || elapsed < t_min);
	t /= n;
	benchmark.print_result(t);
}

// By default long-running timings are disabled (to not annoy the user).
// If the GINAC_RUN_EXPENSIVE_TIMINGS environment variable is set to "1",
// some of them (which are supposed to be relatively fast) will be enabled.
// If GINAC_RUN_EXPENSIVE_TIMINGS is set to "2", all timings are enabled.
static int run_expensive_timings_p()
{
	static int value = 0;
	static int cc = 0;
	static const std::string env_name("GINAC_RUN_EXPENSIVE_TIMINGS");
	if (cc++ == 0) {
		char* envvar = std::getenv(env_name.c_str());
		if (envvar != NULL) {
			value = std::atoi(envvar);
			if (value < 0 || value > 2)
				value = 0;
		}
		if (value) {
			std::cerr << "WARNING: "
				<< "long-running timings are ENABLED."
				<< std::endl
				<< "Unset the \"" << env_name << "\" "
				<< "environment variable skip them."
				<< std::endl;
		}
	}
	return value;
}

#endif // ndef GINAC_CHECK_TEST_RUNNER_H
