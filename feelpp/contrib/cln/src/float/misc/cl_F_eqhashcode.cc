// cl_F equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "base/cl_N.h"
#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/misc/cl_SF_eqhashcode.cc"
#include "float/ffloat/misc/cl_FF_eqhashcode.cc"
#include "float/dfloat/misc/cl_DF_eqhashcode.cc"
#include "float/lfloat/misc/cl_LF_eqhashcode.cc"

namespace cln {

uint32 CL_FLATTEN equal_hashcode (const cl_F& x)
{
	floatcase(x
	,	return equal_hashcode_inline(x);
	,	return equal_hashcode_inline(x);
	,	return equal_hashcode_inline(x);
	,	return equal_hashcode_inline(x);
	);
}

}  // namespace cln
