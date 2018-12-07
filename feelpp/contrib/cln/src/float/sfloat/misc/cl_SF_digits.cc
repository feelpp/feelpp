// float_digits().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

namespace cln {

CL_INLINE uintC CL_INLINE_DECL(float_digits) (const cl_SF& x)
{
	unused x;
	return SF_mant_len+1; // 17
}

}  // namespace cln
