// cl_RA (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_class.h"


// Implementation.

#include "cln/input.h"
#include "cln/rational_io.h"

namespace cln {

cl_read_flags cl_RA_read_flags = {
	syntax_rational,
	lsyntax_all,
	10,
	{ float_format_ffloat, float_format_lfloat_min, true }
};

cl_RA::cl_RA (const char * string)
{
	pointer = as_cl_private_thing(
		read_rational(cl_RA_read_flags,string,NULL,NULL));
}

}  // namespace cln
