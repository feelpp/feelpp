// print_integer().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer_io.h"


// Implementation.

#include "cln/output.h"

namespace cln {

void print_integer (std::ostream& stream, const cl_print_real_flags& flags, const cl_I& z)
{
	print_integer(stream,(const cl_print_rational_flags&)flags,z);
}

}  // namespace cln
