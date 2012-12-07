#include "test_I.h"

int test_I_lognand (int iterations)
{
	int error = 0;
	int i;
	// Check against logand.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(lognand(a,b) == lognot(logand(a,b)), a,b);
	}
	return error;
}
