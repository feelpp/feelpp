#include "test_I.h"

int test_I_abs (int iterations)
{
	int error = 0;
	int i;
	// Check against "-".
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(abs(a) == (a < 0 ? -a : a), a);
	}
	return error;
}
