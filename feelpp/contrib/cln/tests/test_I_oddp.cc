#include "test_I.h"

int test_I_oddp (int iterations)
{
	int error = 0;
	int i;
	// Check against division.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I r = mod(a,2);
		ASSERT1(oddp(a) == (r!=0), a);
	}
	return error;
}
