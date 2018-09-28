#include "test_I.h"

int test_I_logxor (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(logxor(a,b) == logxor(b,a), a,b);
	}
	// Check associativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I c = testrandom_I();
		ASSERT3(logxor(logxor(a,b),c) == logxor(a,logxor(b,c)), a,b,c);
	}
	// Check special cases 0 and -1.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(logxor(a,0) == a, a);
		ASSERT1(logxor(a,a) == 0, a);
		ASSERT1(logxor(a,-1) == lognot(a), a);
		ASSERT1(logxor(a,lognot(a)) == -1, a);
	}
	return error;
}
