// cl_I (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer_class.h"


// Implementation.

#include "cln/input.h"
#include "cln/integer_io.h"

namespace cln {

cl_read_flags cl_I_read_flags = {
	syntax_integer,
	lsyntax_all,
	10,
	{ float_format_ffloat, float_format_lfloat_min, true }
};

cl_I::cl_I (const char * string)
{
	pointer = as_cl_private_thing(
		read_integer(cl_I_read_flags,string,NULL,NULL));
}

}  // namespace cln
