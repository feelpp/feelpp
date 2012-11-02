#include "test_I.h"

int test_I_logcount (int iterations)
{
	int error = 0;
	int i;
	// Check additivity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		if (a >= 0 && b >= 0)
			ASSERT2(logcount(a) == logcount(logand(a,b)) + logcount(logandc2(a,b)), a,b);
	}
	// Check behaviour on sign change.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(logcount(lognot(a)) == logcount(a), a);
	}
	// Check shift invariance.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		sintL b = random32() % 1024;
		ASSERT2(logcount(ash(a,b)) == logcount(a) + (minusp(a) ? b : 0), a,b);
	}
	return error;
}
