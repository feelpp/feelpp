#include "test_MI.h"

int test_MI_mul (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI b = R->canonhom(testrandom_I());
		ASSERT3(a*b == b*a, m,a,b);
	}
	// Check associativity.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI b = R->canonhom(testrandom_I());
		cl_MI c = R->canonhom(testrandom_I());
		ASSERT4((a*b)*c == a*(b*c), m,a,b,c);
	}
	// Check second binomial formula.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI b = R->canonhom(testrandom_I());
		ASSERT3((a+b)*(a-b) == a*a-b*b, m,a,b);
	}
	// Check distributive formula.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI b = R->canonhom(testrandom_I());
		cl_MI c = R->canonhom(testrandom_I());
		ASSERT4((a+c)*(b+c) == a*b+(a+b)*c+c*c, m,a,b,c);
	}
	// Check special cases 0, 1, -1.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		cl_MI z = R->zero();
		cl_MI o = R->one();
		cl_MI mo = R->canonhom(-1);
		ASSERT2(a*z == z, m,a);
		ASSERT2(a*o == a, m,a);
		ASSERT2(a*mo == -a, m,a);
	}
	return error;
}
