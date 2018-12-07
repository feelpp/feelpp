#include "test_I.h"

int test_I_mul (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(a*b == b*a, a,b);
	}
	// Check associativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I c = testrandom_I();
		ASSERT3((a*b)*c == a*(b*c), a,b,c);
	}
	// Check second binomial formula.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2((a+b)*(a-b) == a*a-b*b, a,b);
	}
	// Check distributive formula.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I c = testrandom_I();
		ASSERT3((a+c)*(b+c) == a*b+(a+b)*c+c*c, a,b,c);
	}
	// Check special cases 0, 1, -1.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		ASSERT1(a*0 == 0, a);
		ASSERT1(a*1 == a, a);
		ASSERT1(a*-1 == -a, a);
	}
	return error;
}
