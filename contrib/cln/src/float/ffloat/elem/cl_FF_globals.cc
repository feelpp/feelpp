// Global variables for cl_FF.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat_class.h"
#include "float/ffloat/cl_FF.h"


// Implementation.

namespace cln {

#if !defined(CL_WIDE_POINTERS)

const cl_FF cl_FF_0 = cl_FF_0; // 0.0f0

const cl_FF cl_FF_1 = cl_FF_1; // 1.0f0

const cl_FF cl_FF_minus1 = cl_FF_minus1; // -1.0f0

int cl_FF_globals_init_helper::count = 0;

cl_FF_globals_init_helper::cl_FF_globals_init_helper()
{
	if (count++ == 0) {
		new ((void *)&cl_FF_0) cl_FF(allocate_ffloat(0)); // 0.0f0
		new ((void *)&cl_FF_1) cl_FF(encode_FF(0,1,bit(FF_mant_len))); // 1.0f0
		new ((void *)&cl_FF_minus1) cl_FF(encode_FF(-1,1,bit(FF_mant_len))); // -1.0f0
	}
}

cl_FF_globals_init_helper::~cl_FF_globals_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up
	}
}
#endif

}  // namespace cln

