#include "test_I.h"

int test_I_logtest (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(logtest(a,b) == logtest(b,a), a,b);
	}
	// Check some axioms.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(!logtest(a,logand(b,lognot(a))), a,b);
		ASSERT2(!logtest(logand(b,lognot(a)),a), a,b);
		ASSERT2(!logtest(b,logand(a,lognot(b))), a,b);
		ASSERT2(!logtest(logand(a,lognot(b)),b), a,b);
		if (a != 0)
			ASSERT2(logtest(a,logior(a,b)), a,b);
		if (b != 0)
			ASSERT2(logtest(b,logior(a,b)), a,b);
	}
	// Check some special cases, against ash and logbitp.
	for (i = iterations; i > 0; i--) {
		uintL a = random32() % 1024;
		cl_I b = testrandom_I();
		ASSERT2(logbitp(a,b) == logtest(b,ash(1,(sintL)a)), a,b);
	}
	return error;
}
