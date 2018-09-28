#include "test_I.h"

int test_I_div (int iterations)
{
	int error = 0;
	int i;
	// Check floor.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		if (b != 0) {
			cl_I_div_t qr = floor2(a,b);
			const cl_I& q = qr.quotient;
			const cl_I& r = qr.remainder;
			ASSERT2(a == q*b+r, a,b);
			ASSERT2(b >= 0 ? (r >= 0 && r < b) : (r <= 0 && r > b), a,b);
		}
	}
	// Check ceiling.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		if (b != 0) {
			cl_I_div_t qr = ceiling2(a,b);
			const cl_I& q = qr.quotient;
			const cl_I& r = qr.remainder;
			ASSERT2(a == q*b+r, a,b);
			ASSERT2(b >= 0 ? (r <= 0 && r > -b) : (r >= 0 && r < -b), a,b);
		}
	}
	// Check truncate.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		if (b != 0) {
			cl_I_div_t qr = truncate2(a,b);
			const cl_I& q = qr.quotient;
			const cl_I& r = qr.remainder;
			ASSERT2(a == q*b+r, a,b);
			if (b >= 0)
				if (a >= 0)
					{ ASSERT2(r >= 0 && r < b, a,b); }
				else
					{ ASSERT2(r <= 0 && r > -b, a,b); }
			else
				if (a >= 0)
					{ ASSERT2(r >= 0 && r < -b, a,b); }
				else
					{ ASSERT2(r <= 0 && r > b, a,b); }
		}
	}
	// Check round.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		cl_I b = testrandom_I();
		if (b != 0) {
			cl_I_div_t qr = round2(a,b);
			const cl_I& q = qr.quotient;
			const cl_I& r = qr.remainder;
			ASSERT2(a == q*b+r, a,b);
			ASSERT2(2*abs(r) <= abs(b), a,b);
			if (2*abs(r) == abs(b))
				ASSERT2(evenp(q), a,b);
		}
	}
	return error;
}
