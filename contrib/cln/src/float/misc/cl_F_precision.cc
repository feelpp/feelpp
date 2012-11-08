// float_precision().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline2.h"
#include "float/sfloat/misc/cl_SF_precision.cc"
#include "float/ffloat/misc/cl_FF_precision.cc"
#include "float/dfloat/misc/cl_DF_precision.cc"
#include "float/lfloat/misc/cl_LF_precision.cc"

namespace cln {

uintC CL_FLATTEN float_precision (const cl_F& x)
{
	floatcase(x
	,	return float_precision(x);
	,	return float_precision(x);
	,	return float_precision(x);
	,	return float_precision(x);
	);
}

}  // namespace cln
