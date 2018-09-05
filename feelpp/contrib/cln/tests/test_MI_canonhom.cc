#include "test_MI.h"

int test_MI_canonhom (int iterations)
{
	int error = 0;
	int i;
	// Check canonhom followed by retract.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_I a = testrandom_I();
		ASSERT2(R->retract(R->canonhom(a)) == (m==0 ? a : mod(a,abs(m))), m,a);
	}
	return error;
}
