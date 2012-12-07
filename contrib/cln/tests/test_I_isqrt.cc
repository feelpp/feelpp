#include "test_I.h"

int test_I_isqrt (int iterations)
{
	int error = 0;
	int i;
	// Check against "*".
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		if (a >= 0) {
			cl_I w;
			bool squarep = isqrt(a,&w);
			ASSERT1(w >= 0 && expt_pos(w,2) <= a && a < expt_pos(w+1,2), a);
			ASSERT1(squarep ? w*w==a : w*w<a, a);
		}
	}
	// Check certain special cases.
	for (i = iterations; i > 0; i--) {
		cl_I a = abs(testrandom_I());
		cl_I w;
		// Check a^2 is a square.
		ASSERT1(isqrt(a*a,&w) && w == a, a);
		// Check a^2+1 is not a square, except when a=0.
		ASSERT1((a==0) || (!isqrt(a*a+1,&w) && w == a), a);
		// Check a^2+2*a is not a square, except when a=0.
		ASSERT1(isqrt(a*(a+2),&w)==(a==0) && w == a, a);
	}
	return error;
}
