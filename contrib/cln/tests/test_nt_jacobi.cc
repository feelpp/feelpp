#include "test_nt.h"

int test_nt_jacobi (int iterations)
{
	int error = 0;
	int i;
	// Check special cases.
	for (i = iterations; i > 0; i--) {
		cl_I b = 1+2*abs(testrandom_I());
		cl_I c = testrandom_I();
		if (b==1) {
			ASSERT2(jacobi(b*c,b) == 1, b,c);
		} else {
			ASSERT2(jacobi(b*c,b) == 0, b,c);
		}
	}
	for (i = iterations; i > 0; i--) {
		cl_I b = 1+2*abs(testrandom_I());
		if (mod(b,4)==1) {
			ASSERT1(jacobi(-1,b) == 1, b);
		} else {
			ASSERT1(jacobi(-1,b) == -1, b);
		}
		if (abs(mod(b,8)-4)==1) {
			ASSERT1(jacobi(2,b) == -1, b);
		} else {
			ASSERT1(jacobi(2,b) == 1, b);
		}
	}
	// Check quadratic residues.
	for (i = iterations; i > 0; i--) {
		cl_I b = 1+2*abs(testrandom_I());
		cl_I c = testrandom_I();
		cl_I a = mod(square(c),b);
		if (gcd(a,b)==1) {
			ASSERT2(jacobi(a,b) == 1, b,c);
		} else {
			ASSERT2(jacobi(a,b) == 0, b,c);
		}
	}
	// Check homomorphism (fixed b).
	for (i = iterations; i > 0; i--) {
		cl_I b = 1+2*abs(testrandom_I());
		cl_I a1 = testrandom_I();
		cl_I a2 = testrandom_I();
		ASSERT3(jacobi(a1,b)*jacobi(a2,b) == jacobi(a1*a2,b), b,a1,a2);
	}
	// Check homomorphism (fixed a).
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b1 = 1+2*abs(testrandom_I());
		cl_I b2 = 1+2*abs(testrandom_I());
		ASSERT3(jacobi(a,b1)*jacobi(a,b2) == jacobi(a,b1*b2), a,b1,b2);
	}
	return error;
}
