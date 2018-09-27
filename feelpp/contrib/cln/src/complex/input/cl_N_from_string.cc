// cl_N (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex_class.h"


// Implementation.

#include "cln/input.h"
#include "cln/complex_io.h"

namespace cln {

cl_read_flags cl_N_read_flags = {
	syntax_number,
	lsyntax_all,
	10,
	{ float_format_ffloat, float_format_lfloat_min, true }
};

cl_N::cl_N (const char * string)
{
	pointer = as_cl_private_thing(
		read_complex(cl_N_read_flags,string,NULL,NULL));
}

}  // namespace cln
