// cl_miller_rabin_test().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "numtheory/cl_IF.h"


// Implementation.

#include "cln/modinteger.h"

namespace cln {

bool cl_miller_rabin_test (const cl_I& n, int count, cl_I* factor)
{
	// [Cohen], section 8.2, algorithm 8.2.2.
	var cl_modint_ring R = find_modint_ring(n); // Z/nZ
	var cl_I m = n-1;
	var uintC e = ord2(m);
	m = m>>e;
	// n-1 = 2^e*m
	var cl_MI one = R->one();
	var cl_MI minusone = R->uminus(one);
	for (int i = 0; i < count; i++) {
		// Choosing aa small makes the expt_pos faster.
		var cl_I aa = (i == 0
		               ? (cl_I) 2
		               : i <= cl_small_prime_table_size
		               ? (cl_I) (unsigned int) cl_small_prime_table[i-1] // small prime
		               : 2+random_I(n-2)); // or random >=2, <n
		if (aa >= n)
			break;
		// Now 1 < aa < n.
		var cl_MI a = R->canonhom(aa);
		var cl_MI b = R->expt_pos(a,m); // b = a^m
		if (b == one)
			goto passed;
		for (uintC s = e; s > 0; s--) {
			if (b == minusone)
				goto passed;
			var cl_MI new_b = R->square(b);
			if (new_b == one) {
				// (b-1)*(b+1) == 0 mod n, hence n not prime.
				if (factor)
					*factor = gcd(R->retract(b)-1,n);
				return false;
			}
			b = new_b;
		}
		// b = a^(2^e*m) = a^(n-1), but b != -1 mod n, hence n not prime.
		if (factor) {
			var cl_I g = gcd(aa,n);
			if (g > 1)
				*factor = g;
			else
				*factor = 0;
		}
		return false;
	    passed:
		;
	}
	return true;
}

}  // namespace cln
