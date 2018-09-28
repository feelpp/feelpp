#include "test_MI.h"

int test_MI_expt (int iterations)
{
	int error = 0;
	int i;
	// Check special caSes 0, 1, 2.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_MI a = R->canonhom(testrandom_I());
		ASSERT2(expt(a,0) == R->one(), m,a);
		ASSERT2(expt(a,1) == a, m,a);
		ASSERT2(expt(a,2) == a*a, m,a);
	}
	// Check special cases -1, -2.
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		cl_modint_ring R = find_modint_ring(m);
		cl_I ai = testrandom_I();
		if (gcd(m,ai)==1) {
			cl_MI a = R->canonhom(ai);
			cl_MI ar = R->recip(a);
			ASSERT2(expt(a,-1) == ar, m,a);
			ASSERT2(expt(a,-2) == ar*ar, m,a);
		}
	}
	// Check homomorphism (fixed exponent).
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		if (!zerop(m)) { // avoid generating huge numbers
			cl_modint_ring R = find_modint_ring(m);
			cl_MI a = R->canonhom(testrandom_I());
			cl_MI b = R->canonhom(testrandom_I());
			cl_I e = abs(testrandom_I());
			ASSERT4(expt(a,e)*expt(b,e) == expt(a*b,e), m,a,b,e);
		}
	}
	// Check distributive formulas (fixed base).
	for (i = iterations; i > 0; i--) {
		cl_I m = testrandom_I();
		if (!zerop(m)) { // avoid generating huge numbers
			cl_modint_ring R = find_modint_ring(m);
			cl_MI a = R->canonhom(testrandom_I());
			cl_I e = abs(testrandom_I());
			cl_I f = abs(testrandom_I());
			ASSERT4(expt(a,e)*expt(a,f) == expt(a,e+f), m,a,e,f);
			ASSERT4(expt(expt(a,e),f) == expt(a,e*f), m,a,e,f);
		}
	}
	return error;
}
