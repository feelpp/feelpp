// float_precision().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

#include "base/cl_inline.h"
#include "float/ffloat/elem/cl_FF_zerop.cc"

namespace cln {

CL_INLINE2 uintC CL_INLINE2_DECL(float_precision) (const cl_FF& x)
{
	if (zerop_inline(x)) return 0;
	return FF_mant_len+1; // 24
}

}  // namespace cln
