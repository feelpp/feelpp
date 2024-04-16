#include "test_I.h"
#include <cln/input.h>
#include <sstream>

int test_I_io (int iterations)
{
	int error = 0;
	for (int i = iterations; i > 0; i--) {
		cl_I a = testrandom_I();
		unsigned base = iterations % (36-1) + 2;
		cl_read_flags rflags = {syntax_integer, lsyntax_standard, base};
		stringstream buf;
		print_integer(buf, base, a);
		cl_I b;
		try {
			b = read_integer(buf, rflags);
			ASSERT1(a == b, a);
		} catch (runtime_exception& err) {
			std::cerr << "Got an error while parsing \"" 
				<< buf.str() << "\" with base = " << base
				<< " (in decimal: " << a << ")" << std::endl;
			std::cerr << "Details: " << err.what() << std::endl;
			++error;
			break;
		}

	}
	return error;
}
