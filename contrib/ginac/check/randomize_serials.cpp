/** @file randomize_serials.cpp
 *
 *  Utilitiy function used by the benchmarks.
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

#include "ginac.h"
using namespace GiNaC;

#include <cstdlib>
#include <ctime>
using namespace std;

/** Generate a random amount of symbols and destroy them again immediatly.
 *  This operation effectively makes the serial numbers of all subsequent
 *  symbols unpredictable.  If the serials are unpredictable, then so are
 *  their hash values.  If the hash values are unpredictable, then so are
 *  the canonical orderings.  If the canonical orderings are unpredictable,
 *  all subsequent times are subject to some variation.  This variation,
 *  however is natural and desireable for two reasons: First, we cannot know
 *  how many symbols have been generated before in real world computations.
 *  Second, the following timings are subject to some semi-random variation
 *  anyways because short timings need to be repeated until enough time has
 *  gone by for the measurement to be reliable.  During this process the serial
 *  numbers will be shifted anyways in a semi-random way.  It is better not
 *  to lull the user in a false sense of reproducibility and instead confront
 *  her with the normal variation to be expected.
 */
void randomify_symbol_serials()
{
	srand(time(NULL));
	const int m = rand() % 666;
	for (int s=0; s<m; ++s ) {
		symbol("dummy");
	}
}
