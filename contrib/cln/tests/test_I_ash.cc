#include "test_I.h"

int test_I_ash (int iterations)
{
	int error = 0;
	int i;
	// Check against "*" and floor1.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		uintL b = random32() % 1024;
		cl_I pow2 = 1;
		for (uintL j = 0; j < b; j++)
			pow2 = pow2+pow2;
		ASSERT1(ash(1,(sintL)b) == pow2, b);
		ASSERT1(ash(1,(cl_I)b) == pow2, b);
		ASSERT2(ash(a,(sintL)b) == a*pow2, a,b);
		ASSERT2(ash(a,-(sintL)b) == floor1(a,pow2), a,b);
	}
	// Check homomorphism w.r.t. second argument.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		sintL b = random32() % 1024;
		sintL c = random32() % 1024;
		ASSERT3(ash(ash(a,b),c) == ash(a,b+c), a,b,c);
		ASSERT3(ash(ash(a,b),-c) == ash(a,b-c), a,b,c);
		ASSERT3(ash(ash(a,-b),-c) == ash(a,-b-c), a,b,c);
	}
	// Check against each other.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		uintL b = random32() % 1024;
		ASSERT2(ash(a,(cl_I)b) == ash(a,(sintL)b), a,b);
		ASSERT2(ash(a,-(cl_I)b) == ash(a,-(sintL)b), a,b);
	}
	return error;
}
