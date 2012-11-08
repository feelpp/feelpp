// plusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline2.h"
#include "float/sfloat/elem/cl_SF_plusp.cc"
#include "float/ffloat/elem/cl_FF_plusp.cc"
#include "float/dfloat/elem/cl_DF_plusp.cc"
#include "float/lfloat/elem/cl_LF_plusp.cc"

namespace cln {

bool CL_FLATTEN plusp (const cl_F& x)
{
	floatcase(x
	,	return plusp_inline(x);
	,	return plusp_inline(x);
	,	return plusp_inline(x);
	,	return plusp_inline(x);
	);
}

}  // namespace cln
