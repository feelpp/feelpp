#include "test_I.h"

int test_I_minus (int iterations)
{
	int error = 0;
	int i;
	// Check anti-commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2((a-b) + (b-a) == 0, a,b);
	}
	// Check associativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I c = testrandom_I();
		ASSERT3(a-(b-c) == c-(b-a), a,b,c);
	}
	// Check special case 0.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(a-0 == a, a);
		ASSERT1(0-(0-a) == a, a);
	}
	return error;
}
