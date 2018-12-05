#include "test_I.h"

int test_I_lognor (int iterations)
{
	int error = 0;
	int i;
	// Check against logior.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(lognor(a,b) == lognot(logior(a,b)), a,b);
	}
	return error;
}
