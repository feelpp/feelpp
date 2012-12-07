// cl_F (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float_class.h"


// Implementation.

#include "cln/float.h"
#include "cln/input.h"
#include "cln/float_io.h"

namespace cln {

cl_read_flags cl_F_read_flags = {
	syntax_float,
	lsyntax_all,
	10,
	{ float_format_ffloat, float_format_lfloat_min, true }
};

cl_F::cl_F (const char * string)
{
	pointer = as_cl_private_thing(
		read_float(cl_F_read_flags,string,NULL,NULL));
}

}  // namespace cln
