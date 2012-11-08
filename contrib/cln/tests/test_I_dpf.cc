#include "test_I.h"

int test_I_deposit_field (int iterations)
{
	int error = 0;
	int i;
	// Check against ash.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		sintL s = random32() % 1024;
		sintL p = random32() % 1024;
		cl_I mask = ash(ash(1,s)-1,p);
		ASSERT4(deposit_field(a,b,cl_byte(s,p)) == logxor(logand(a,mask),logandc2(b,mask)), a,s,p,b);
	}
	return error;
}
