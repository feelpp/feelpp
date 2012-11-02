// Global variables for cl_DF.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/dfloat/cl_DF.h"


// Implementation.

namespace cln {

const cl_DF cl_DF_0 = cl_DF_0;
const cl_DF cl_DF_1 = cl_DF_1;
const cl_DF cl_DF_minus1 = cl_DF_minus1;

int cl_DF_globals_init_helper::count = 0;

cl_DF_globals_init_helper::cl_DF_globals_init_helper()
{
	if (count++ == 0) {
#if (cl_word_size == 64)
		new ((void *)&cl_DF_0) cl_DF(allocate_dfloat(0)); // 0.0d0
		new ((void *)&cl_DF_1) cl_DF(encode_DF(0, 1, bit(DF_mant_len))); // 1.0d0
		new ((void *)&cl_DF_minus1) cl_DF(encode_DF(-1,1,bit(DF_mant_len))); // -1.0d0
#else
		new ((void *)&cl_DF_0) cl_DF(allocate_dfloat(0, 0)); // 0.0d0
		new ((void *)&cl_DF_1) cl_DF(encode_DF(0, 1, bit(DF_mant_len - 32), 0)); // 1.0d0
		new ((void *)&cl_DF_minus1) cl_DF(encode_DF(-1, 1, bit(DF_mant_len - 32), 0)); // -1.0d0
#endif
	}
}
cl_DF_globals_init_helper::~cl_DF_globals_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up
	}
}

}  // namespace cln

