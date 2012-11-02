// float_digits().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

CL_INLINE uintC CL_INLINE_DECL(float_digits) (const cl_DF& x)
{
	unused x;
	return DF_mant_len+1; // 53
}

}  // namespace cln
