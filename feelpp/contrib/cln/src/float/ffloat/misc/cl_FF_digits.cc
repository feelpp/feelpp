// float_digits().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

CL_INLINE uintC CL_INLINE_DECL(float_digits) (const cl_FF& x)
{
	unused x;
	return FF_mant_len+1; // 24
}

}  // namespace cln
