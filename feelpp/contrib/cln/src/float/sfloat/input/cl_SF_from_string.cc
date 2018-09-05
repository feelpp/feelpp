// cl_SF (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat_class.h"


// Implementation.

#include "cln/sfloat.h"
#include "cln/input.h"
#include "cln/float_io.h"

namespace cln {

cl_read_flags cl_SF_read_flags = {
	syntax_sfloat,
	lsyntax_all,
	10,
	{ float_format_sfloat, float_format_lfloat_min, false }
};

cl_SF::cl_SF (const char * string)
{
	pointer = as_cl_private_thing(
		As(cl_SF)(read_float(cl_SF_read_flags,string,NULL,NULL)));
}

}  // namespace cln
