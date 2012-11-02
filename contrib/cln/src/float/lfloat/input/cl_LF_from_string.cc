// cl_LF (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat_class.h"


// Implementation.

#include "cln/lfloat.h"
#include "cln/input.h"
#include "cln/float_io.h"

namespace cln {

cl_read_flags cl_LF_read_flags = {
	syntax_lfloat,
	lsyntax_all,
	10,
	{ float_format_lfloat_min, float_format_lfloat_min, false }
};

cl_LF::cl_LF (const char * string)
{
	pointer = as_cl_private_thing(
		As(cl_LF)(read_float(cl_LF_read_flags,string,NULL,NULL)));
}

}  // namespace cln
