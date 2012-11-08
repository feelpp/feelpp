#include "test_I.h"

int test_I_integer_length (int iterations)
{
	int error = 0;
	int i;
	// Check against ash.
	for (i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		uintC l = integer_length(a);
		if (a >= 0) {
			int b = 0;
			if (a < ash(1,l)) {
				if (a == 0)
					b = (l == 0);
				else
					b = (l > 0 && a >= ash(1,l-1));
			}
			ASSERT1(b, a);
		} else {
			int b = 0;
			if (a >= ash(-1,l)) {
				if (a == -1)
					b = (l == 0);
				else
					b = (l > 0 && a < ash(-1,l-1));
			}
			ASSERT1(b, a);
		}
	}
	return error;
}
