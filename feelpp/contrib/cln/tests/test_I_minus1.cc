#include "test_I.h"

int test_I_minus1 (int iterations)
{
	int error = 0;
	int i;
	// Check against "+".
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(minus1(a) == -1+a, a);
	}
	return error;
}
