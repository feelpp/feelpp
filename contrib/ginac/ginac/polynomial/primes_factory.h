/** @file primes_factory.h
 *
 *  Factory for prime numbers. */

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

#ifndef GINAC_CHINREM_GCD_PRIMES_FACTORY_H
#define GINAC_CHINREM_GCD_PRIMES_FACTORY_H

#include "smod_helpers.h"
#include "debug.h"

#include <cln/integer.h>
#include <cln/numtheory.h>
#include <limits>

namespace GiNaC {

/**
 * Find a `big' prime p such that lc mod p != 0. Helper class used by modular
 * GCD algorithm.
 */
class primes_factory
{
private:
	// These primes need to be large enough, so that the number of images
	// we need to reconstruct the GCD (in Z) is reasonable. On the other
	// hand, they should be as small as possible, so that operations on
	// coefficients are efficient. Practically this means we coefficients
	// should be native integers. (N.B.: as of now chinrem_gcd uses cl_I
	// or even numeric. Eventually this will be fixed).
	cln::cl_I last;
	// This ensures coefficients are immediate.
	static const int immediate_bits = 8*sizeof(void *) - __alignof__(void *);
	static const long opt_hint = (1L << (immediate_bits >> 1)) - 1;
public:
	primes_factory()
	{
		last = cln::nextprobprime(cln::cl_I(opt_hint));
	}

	bool operator()(long& p, const cln::cl_I& lc)
	{
		static const cln::cl_I maxval(std::numeric_limits<long>::max());
		while (last < maxval) {
			long p_ = cln::cl_I_to_long(last);
			last = cln::nextprobprime(last + 1);

			if (!zerop(smod(lc, p_))) {
				p = p_;
				return true;
			}
		}
		return false;
	}

	bool has_primes() const
	{
		static const cln::cl_I maxval(std::numeric_limits<long>::max());
		return last < maxval;
	}
};

} // namespace GiNaC

#endif // ndef GINAC_CHINREM_GCD_PRIMES_FACTORY_H
