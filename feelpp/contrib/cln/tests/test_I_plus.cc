#include "test_I.h"

int test_I_plus (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(a+b == b+a, a,b);
	}
	// Check associativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I c = testrandom_I();
		ASSERT3((a+b)+c == a+(b+c), a,b,c);
	}
	// Check special case 0.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(a+0 == a, a);
	}
	return error;
}
