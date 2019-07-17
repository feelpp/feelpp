// float_precision().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

#include "base/cl_inline.h"
#include "float/dfloat/elem/cl_DF_zerop.cc"

namespace cln {

CL_INLINE2 uintC CL_INLINE2_DECL(float_precision) (const cl_DF& x)
{
	if (zerop_inline(x)) return 0;
	return DF_mant_len+1; // 53
}

}  // namespace cln
