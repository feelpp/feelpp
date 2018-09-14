// Global variables for cl_LF.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/number.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"

namespace cln {

// Only needed for the default constructor of cl_LF.
const cl_LF cl_LF_0 = cl_LF_0; // 0.0L0

int cl_LF_globals_init_helper::count = 0;
cl_LF_globals_init_helper::cl_LF_globals_init_helper()
{
	if (count++ == 0)
		new ((void *)&cl_LF_0) cl_LF(encode_LF0(LF_minlen)); // 0.0L0
}

cl_LF_globals_init_helper::~cl_LF_globals_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up
	}
}

}  // namespace cln

