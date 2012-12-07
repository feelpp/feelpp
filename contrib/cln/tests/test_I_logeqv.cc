#include "test_I.h"

int test_I_logeqv (int iterations)
{
	int error = 0;
	int i;
	// Check against logxor.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(logeqv(a,b) == lognot(logxor(a,b)), a,b);
	}
	return error;
}
