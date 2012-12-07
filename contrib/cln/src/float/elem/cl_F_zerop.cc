// zerop().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_zerop.cc"
#include "float/ffloat/elem/cl_FF_zerop.cc"
#include "float/dfloat/elem/cl_DF_zerop.cc"
#include "float/lfloat/elem/cl_LF_zerop.cc"

namespace cln {

bool CL_FLATTEN zerop (const cl_F& x)
{
	floatcase(x
	,	return zerop_inline(x);
	,	return zerop_inline(x);
	,	return zerop_inline(x);
	,	return zerop_inline(x);
	);
}

}  // namespace cln
