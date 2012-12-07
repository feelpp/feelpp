// cl_trialdivision().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "numtheory/cl_IF.h"


// Implementation.

#include "base/cl_low.h"

namespace cln {

uint32 cl_trialdivision (uint32 n, uint32 d1, uint32 d2)
{
	var uintL i = cl_small_prime_table_search(d1);
	var const uint16 * ptr = &cl_small_prime_table[i];
	var const uint16 * ptr_limit = &cl_small_prime_table[cl_small_prime_table_search(d2+1)];
	for ( ; ptr < ptr_limit; ptr++) {
		var uint32 prime = *ptr;
		var uint32 r;
		r = n % prime; // or: divu_3232_3232(n,prime,,r=);
		if (r == 0)
			return prime;
	}
	return 0;
}

}  // namespace cln
