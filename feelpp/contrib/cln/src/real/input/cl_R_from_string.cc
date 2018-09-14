// cl_R (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real_class.h"


// Implementation.

#include "cln/input.h"
#include "cln/real_io.h"

namespace cln {

cl_read_flags cl_R_read_flags = {
	syntax_real,
	lsyntax_all,
	10,
	{ float_format_ffloat, float_format_lfloat_min, true }
};

cl_R::cl_R (const char * string)
{
	pointer = as_cl_private_thing(
		read_real(cl_R_read_flags,string,NULL,NULL));
}

}  // namespace cln
