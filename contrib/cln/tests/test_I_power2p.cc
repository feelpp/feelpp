#include "test_I.h"

int test_I_power2p (int iterations)
{
	int error = 0;
	int i;
	// Check powers of 2.
	for (i = iterations; i > 0; i--) {
		uintL n = random32() % 1024;
		cl_I a = ash(1,n);
		ASSERT1(power2p(a) == n+1, n);
	}
	// Check against logcount.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		if (a > 0)
			ASSERT1(power2p(a) == (logcount(a) == 1 ? integer_length(a) : 0), a);
	}
	return error;
}
