#include "test_MI.h"

int test_MI_minus (int iterations)
{
	int error = 0;
	int i;
	// Check anti-commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI b = R->canonhom(testrandom_I());
		cl_MI z = R->zero();
		ASSERT3((a-b) + (b-a) == z, m,a,b);
	}
	// Check associativity.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI b = R->canonhom(testrandom_I());
		cl_MI c = R->canonhom(testrandom_I());
		ASSERT4(a-(b-c) == c-(b-a), m,a,b,c);
	}
	// Check special case 0.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI z = R->zero();
		ASSERT2(a-z == a, m,a);
		ASSERT2(z-(z-a) == a, m,a);
	}
	return error;
}
