#include "test_I.h"

int test_I_lognot (int iterations)
{
	int error = 0;
	int i;
	// Check involution, sign, and against "+".
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = lognot(a);
		ASSERT1(lognot(b) == a, a);
		ASSERT1(minusp(a) != minusp(b), a);
		ASSERT1(a+b == -1, a);
	}
	// Check special cases 0 and -1.
	ASSERT(lognot(0) == -1);
	return error;
}
