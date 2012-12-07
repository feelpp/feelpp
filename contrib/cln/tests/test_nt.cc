#include <iostream>

extern int test_nt_jacobi (int iterations);

#define RUN(tester,iterations)  \
	std::cout << "Testing "#tester"..." << std::endl; \
	error |= tester (iterations);

int test_nt (int iterations)
{
	int error = 0;
	RUN(test_nt_jacobi,iterations);
	return error;
}
