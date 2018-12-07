// float_exponent().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/misc/cl_SF_exponent.cc"
#include "float/ffloat/misc/cl_FF_exponent.cc"
#include "float/dfloat/misc/cl_DF_exponent.cc"
#include "float/lfloat/misc/cl_LF_exponent.cc"

namespace cln {

sintE CL_FLATTEN float_exponent (const cl_F& x)
{
	floatcase(x
	,	return float_exponent_inline(x);
	,	return float_exponent_inline(x);
	,	return float_exponent_inline(x);
	,	return float_exponent_inline(x);
	);
}

}  // namespace cln
