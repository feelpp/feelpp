// cl_FF (const char *) constructor.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat_class.h"


// Implementation.

#include "cln/ffloat.h"
#include "cln/input.h"
#include "cln/float_io.h"

namespace cln {

cl_read_flags cl_FF_read_flags = {
	syntax_ffloat,
	lsyntax_all,
	10,
	{ float_format_ffloat, float_format_lfloat_min, false }
};

cl_FF::cl_FF (const char * string)
{
	pointer = as_cl_private_thing(
		As(cl_FF)(read_float(cl_FF_read_flags,string,NULL,NULL)));
}

}  // namespace cln
