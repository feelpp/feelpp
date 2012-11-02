#include <iostream>
#include "base/cl_macros.h"

extern int test_MI_canonhom (int iterations);
extern int test_MI_plus (int iterations);
extern int test_MI_minus (int iterations);
extern int test_MI_mul (int iterations);
extern int test_MI_recip (int iterations);
extern int test_MI_div (int iterations);
extern int test_MI_expt (int iterations);

#define RUN(tester,iterations)  \
	std::cout << "Testing "#tester"..." << std::endl; \
	error |= tester (iterations);

int test_MI (int iterations)
{
	int error = 0;
	RUN(test_MI_canonhom,iterations);
	RUN(test_MI_plus,iterations);
	RUN(test_MI_minus,iterations);
	RUN(test_MI_mul,iterations);
	RUN(test_MI_recip,iterations);
	RUN(test_MI_div,iterations);
	RUN(test_MI_expt,ceiling(iterations,20));
	return error;
}
