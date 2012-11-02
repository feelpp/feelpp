#include "test_I.h"
#define floor(a,b)  ((a) / (b))

int test_I_sqrtp (int iterations)
{
	int error = 0;
	int i;
	// Check against isqrt.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		if (a >= 0) {
			cl_I w;
			bool squarep = sqrtp(a,&w);
			cl_I correct_w;
			bool correct_squarep = isqrt(a,&correct_w);
			ASSERT1(squarep == correct_squarep, a);
			if (squarep)
				ASSERT1(w == correct_w, a);
		}
	}
	// Check certain special cases.
	for (i = iterations; i > 0; i--) {
		cl_I a = abs(testrandom_I());
		cl_I w;
		// Check a^2 is a square.
		ASSERT1(sqrtp(a*a,&w) && w == a, a);
		// Check a^2+1 is not a square, except when a=0.
		if (a > 0) ASSERT1(!isqrt(a*a+1,&w) && w == a, a);
		// Check a^2+2*a is not a square, except when a=0.
		ASSERT1(isqrt(a*(a+2),&w)==(a==0) && w == a, a);
	}
	return error;
}
