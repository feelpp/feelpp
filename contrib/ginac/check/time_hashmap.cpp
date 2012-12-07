/** @file time_hashmap.cpp
 *
 *  Timings for exhashmap<> operations. */

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

template <class T>
static void run_timing(unsigned size, double &time_insert, double &time_find, double &time_erase)
{
	vector<symbol> S;
	T M;
	timer t;

	S.reserve(size);
	for (unsigned i=0; i<size; ++i)
		S.push_back(symbol());

	t.start();
	for (unsigned i=0; i<size; ++i)
		M[S[i]] = S[(i+1)%size];
	time_insert = t.read();

	t.start();
	for (unsigned i=0; i<size; ++i) {
		if (!M[S[i]].is_equal(S[(i+1)%size])) {
			clog << "map lookup failed" << endl;
			return;
		}
	}
	time_find = t.read();

	t.start();
	for (unsigned i=0; i<size; ++i) {
		if (M.erase(S[i]) != 1) {
			clog << "erasing element " << S[i] << " failed" << endl;
			return;
		}
	}
	if (!M.empty()) {
		clog << "map not empty (size = " << M.size() << ") after erasing all elements" << endl;
		return;
	}
	time_erase = t.read();
}


unsigned time_hashmap()
{
	unsigned result = 0;

	cout << "timing hash map operations" << flush;

	unsigned s[] = {10000, 50000, 100000, 500000};
	vector<unsigned> sizes(s, s+sizeof(s)/sizeof(*s));

	vector<double> times_insert, times_find, times_erase;

	for (vector<unsigned>::const_iterator i = sizes.begin(); i != sizes.end(); ++i) {
		double time_insert, time_find, time_erase;

		run_timing< exhashmap<ex> >(*i, time_insert, time_find, time_erase);

// If you like, you can compare it with this:
//		run_timing< std::map<ex, ex, ex_is_less> >(*i, time_insert, time_find, time_erase);

		times_insert.push_back(time_insert);
		times_find.push_back(time_find);
		times_erase.push_back(time_erase);
		cout << '.' << flush;
	}

	// print the report:
	cout << endl << "          size:\t";
	copy(sizes.begin(), sizes.end(), ostream_iterator<unsigned>(cout, "\t"));
	cout << endl << "      insert/s:\t";
	copy(times_insert.begin(), times_insert.end(), ostream_iterator<double>(cout, "\t"));
	cout << endl << "        find/s:\t";
	copy(times_find.begin(), times_find.end(), ostream_iterator<double>(cout, "\t"));
	cout << endl << "       erase/s:\t";
	copy(times_erase.begin(), times_erase.end(), ostream_iterator<double>(cout, "\t"));
	cout << endl;

	return result;
}

extern void randomify_symbol_serials();

int main(int argc, char** argv)
{
	randomify_symbol_serials();
	cout << setprecision(2) << showpoint;
	return time_hashmap();
}
