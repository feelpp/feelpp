#include "test_I.h"

int test_I_gcd (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(gcd(a,b) == gcd(b,a), a,b);
	}
	// Check some axioms.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I g = gcd(a,b);
		if (g > 1) {
			ASSERT2(mod(a,g) == 0, a,b);
			ASSERT2(mod(b,g) == 0, a,b);
			ASSERT2(gcd(exquo(a,g),exquo(b,g)) == 1, a,b);
		}
		cl_I c = testrandom_I();
		ASSERT3(gcd(a+b*c,b) == g, a,b,c);
	}
	return error;
}
