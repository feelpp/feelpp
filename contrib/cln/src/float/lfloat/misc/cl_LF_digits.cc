// float_digits().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"

namespace cln {

CL_INLINE uintC CL_INLINE_DECL(float_digits) (const cl_LF& x)
{
	return intDsize*(uintC)(TheLfloat(x)->len);
}

}  // namespace cln
