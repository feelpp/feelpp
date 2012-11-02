// float_precision().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_zerop.cc"

namespace cln {

CL_INLINE2 uintC CL_INLINE2_DECL(float_precision) (const cl_SF& x)
{
	if (zerop_inline(x)) return 0;
	return SF_mant_len+1; // 17
}

}  // namespace cln
