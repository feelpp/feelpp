// float_exponent().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"

namespace cln {

CL_INLINE sintE CL_INLINE_DECL(float_exponent) (const cl_LF& x)
{
	var uintE uexp = TheLfloat(x)->expo;
	if (uexp==0) { return 0; }
	return (sintE)(uexp - LF_exp_mid);
}

}  // namespace cln
