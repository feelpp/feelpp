// print_rational().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_io.h"


// Implementation.

#include "cln/output.h"

namespace cln {

void print_rational (std::ostream& stream, const cl_print_flags& flags, const cl_RA& z)
{
	print_rational(stream,(const cl_print_number_flags&)flags,z);
}

}  // namespace cln
