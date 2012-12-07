#include "test_I.h"

int test_I_ldb (int iterations)
{
	int error = 0;
	int i;
	// Check against ash.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		sintL s = random32() % 1024;
		sintL p = random32() % 1024;
		ASSERT3(ldb(a,cl_byte(s,p)) == logand(ash(a,-p),ash(1,s)-1), a,s,p);
	}
	return error;
}
