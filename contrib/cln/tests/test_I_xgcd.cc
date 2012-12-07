#include "test_I.h"

int test_I_xgcd (int iterations)
{
	int error = 0;
	int i;
	// Check against gcd.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		cl_I u, v;
		cl_I g = xgcd(a,b,&u,&v);
		ASSERT3(g == gcd(a,b), a,b,g);
		ASSERT4(g == u*a+v*b, a,b,u,v);
		if (a != 0 && b != 0) {
			if (abs(a) == abs(b)) {
				ASSERT4((u == signum(a) && v == 0) || (u == 0 && v == signum(b)), a,b,u,v);
			}
			else if (mod(abs(a),abs(b)) == 0) {
				ASSERT4(u == 0 && v == signum(b), a,b,u,v);
			}
			else if (mod(abs(b),abs(a)) == 0) {
				ASSERT4(u == signum(a) && v == 0, a,b,u,v);
			}
			else {
				ASSERT4(abs(u) <= floor1(abs(b),2*g) && abs(v) <= floor1(abs(a),2*g), a,b,u,v);
			}
		}
	}
	return error;
}
