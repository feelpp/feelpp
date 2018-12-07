// print_real().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real_io.h"


// Implementation.

#include "cln/output.h"

namespace cln {

void print_real (std::ostream& stream, const cl_print_flags& flags, const cl_R& z)
{
	print_real(stream,(const cl_print_number_flags&)flags,z);
}

}  // namespace cln
