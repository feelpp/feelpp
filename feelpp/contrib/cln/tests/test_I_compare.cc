#include "test_I.h"

int test_I_compare (int iterations)
{
	int error = 0;
	int i;
	// Check anticommutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(compare(a,b) == -compare(b,a), a,b);
	}
	// Check  a < b  <==>  a+c < b+c .
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I c = testrandom_I();
		ASSERT3(compare(a,b) == compare(a+c,b+c), a,b,c);
	}
	return error;
}
