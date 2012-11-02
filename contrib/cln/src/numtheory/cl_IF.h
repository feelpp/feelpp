// Integer factorization and primality testing.

#ifndef _CL_IF_H
#define _CL_IF_H

#include "cln/number.h"
#include "cln/integer.h"

namespace cln {

// Table of primes > 2, < 2^16
const uint32 cl_small_prime_table_limit = 65536;
const int cl_small_prime_table_size = 6541;
extern uint16 cl_small_prime_table[cl_small_prime_table_size];

// Given 0 < d <= cl_small_prime_table_limit, return the smallest index i
// such that  cl_small_prime_table[i] >= d.  (Or i = cl_small_prime_table_size
// if none exists.)
inline uintL cl_small_prime_table_search (uint32 d)
{
	var uintL i1 = 0;
	var uintL i2 = cl_small_prime_table_size;
	if (cl_small_prime_table[i1] >= d)
		return i1;
	loop {
		// Here i1 < i2 and
		// cl_small_prime_table[i1] < d <= cl_small_prime_table[i2].
		var uintL i3 = floor(i1+i2,2);
		if (i3 == i1) // (i2-i1 == 1) ?
			return i2;
		if (cl_small_prime_table[i3] >= d)
			i2 = i3;
		else
			i1 = i3;
	}
}

// Trial division.
// Divides n > 0 by the primes in the range d1 <= d <= d2
// (0 < d1 <= d2 <= min(isqrt(n),cl_small_prime_table_limit))
// and returns the divisor d if found, or 0 if no divisor found.
extern uint32 cl_trialdivision (uint32 n, uint32 d1, uint32 d2);
extern uint32 cl_trialdivision (uint32 nhi, uint32 nlo, uint32 d1, uint32 d2);
extern uint32 cl_trialdivision (const cl_I& n, uint32 d1, uint32 d2);

// Miller-Rabin compositeness test.
// Performs count times the Miller-Rabin test on n > 1 odd.
// Returns true if n looks like a prime (with error probability < 4^-count).
// Returns false if n is definitely composite, and then sets factor = some
// nontrivial factor or 0.
extern bool cl_miller_rabin_test (const cl_I& n, int count, cl_I* factor);

}  // namespace cln

#endif /* _CL_IF_H */
