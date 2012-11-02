#include "test_I.h"

int test_I_plus1 (int iterations)
{
	int error = 0;
	int i;
	// Check against "+".
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(plus1(a) == 1+a, a);
	}
	return error;
}
