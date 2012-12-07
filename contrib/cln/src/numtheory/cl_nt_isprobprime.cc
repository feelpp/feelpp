// isprobprime().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/numtheory.h"


// Implementation.

#include "numtheory/cl_IF.h"
#include "cln/integer_io.h"
#include "cln/exception.h"
#include <sstream>

namespace cln {

bool isprobprime (const cl_I& n)
{
	if (!(n > 0)) {
		std::ostringstream buf;
		fprint(buf, n);
		fprint(buf, " is not a positive integer.");
		throw runtime_exception(buf.str());
	}
	// With a Miller-Rabin count = 50 the final error probability is
	// 4^-50 < 10^-30.
	var int count = 50;
	// Step 1: Trial division (rules out 87% of all numbers quickly).
	const uint32 trialdivide_limit = 70;
	var uintC l = integer_length(n);
	if (l <= 32) {
		var uint32 nn = cl_I_to_UL(n);
		if (nn <= cl_small_prime_table_limit) {
			// Table lookup.
			var uintL i = cl_small_prime_table_search(nn);
			if (i < cl_small_prime_table_size
			    && ((unsigned int) cl_small_prime_table[i] == nn
			        || nn == 2))
				return true;
			else
				return false;
		}
		if ((nn % 2) == 0 || cl_trialdivision(nn,1,trialdivide_limit))
			return false;
		// For small n, only few Miller-Rabin tests are needed.
		if (nn < 2000U) count = 1; // {2}
		else if (nn < 1300000U) count = 2; // {2,3}
		else if (nn < 25000000U) count = 3; // {2,3,5}
		else if (nn < 3200000000U) count = 4; // {2,3,5,7}
	} else if (l <= 64) {
		var uint32 nhi = cl_I_to_UL(ldb(n,cl_byte(32,32)));
		var uint32 nlo = cl_I_to_UL(ldb(n,cl_byte(32,0)));
		if ((nlo % 2) == 0 || cl_trialdivision(nhi,nlo,1,trialdivide_limit))
			return false;
	} else {
		if (evenp(n) || cl_trialdivision(n,1,trialdivide_limit))
			return false;
	}
	// Step 2: Miller-Rabin test.
	return cl_miller_rabin_test(n,count,NULL);
}

}  // namespace cln
