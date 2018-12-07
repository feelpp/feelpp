// minusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_minusp.cc"
#include "float/ffloat/elem/cl_FF_minusp.cc"
#include "float/dfloat/elem/cl_DF_minusp.cc"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

bool CL_FLATTEN minusp (const cl_F& x)
{
	floatcase(x
	,	return minusp_inline(x);
	,	return minusp_inline(x);
	,	return minusp_inline(x);
	,	return minusp_inline(x);
	);
}

}  // namespace cln
