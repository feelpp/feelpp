#include "test_MI.h"

int test_MI_recip (int iterations)
{
	int error = 0;
	int i;
	// Check against multiplication.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_I ai = testrandom_I();
		if (gcd(m,ai)==1) {
			cl_MI a = R->canonhom(ai);
			ASSERT2(a*R->recip(a) == R->one(), m,a);
		}
	}
	return error;
}
