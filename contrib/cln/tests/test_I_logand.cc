#include "test_I.h"

int test_I_logand (int iterations)
{
	int error = 0;
	int i;
	// Check commutativity.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(logand(a,b) == logand(b,a), a,b);
	}
	// Check against ash and oddp.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		ASSERT2(oddp(logand(a,b)) == (oddp(a) & oddp(b)), a,b);
		sintL c = random32() % 1024;
		ASSERT3(logand(ash(a,c),ash(b,c)) == ash(logand(a,b),c), a,b,c);
		ASSERT3(logand(ash(a,-c),ash(b,-c)) == ash(logand(a,b),-c), a,b,c);
	}
	return error;
}
