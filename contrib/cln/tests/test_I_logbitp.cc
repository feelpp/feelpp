#include "test_I.h"

int test_I_logbitp (int iterations)
{
	int error = 0;
	int i;
	// Check against ash and oddp.
	for (i = iterations; i > 0; i--) {
		uintL a = random32() % 1024;
		cl_I b = testrandom_I();
		ASSERT2(logbitp(a,b) == oddp(ash(b,-(sintL)a)), a,b);
	}
	// Check against ash and logand.
	for (i = iterations; i > 0; i--) {
		uintL a = random32() % 1024;
		cl_I b = testrandom_I();
		ASSERT2(logbitp(a,b) == !zerop(logand(b,ash(1,(sintL)a))), a,b);
	}
	// Check against each other.
	for (i = iterations; i > 0; i--) {
		uintL a = random32() % 1024;
		cl_I b = testrandom_I();
		ASSERT2(logbitp((cl_I)a,b) == logbitp(a,b), a,b);
	}
	return error;
}
