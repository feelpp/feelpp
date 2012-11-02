// cl_DF (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat_class.h"


// Implementation.

#include "cln/dfloat.h"
#include "cln/input.h"
#include "cln/float_io.h"

namespace cln {

cl_read_flags cl_DF_read_flags = {
	syntax_dfloat,
	lsyntax_all,
	10,
	{ float_format_dfloat, float_format_lfloat_min, false }
};

cl_DF::cl_DF (const char * string)
{
	pointer = as_cl_private_thing(
		As(cl_DF)(read_float(cl_DF_read_flags,string,NULL,NULL)));
}

}  // namespace cln
