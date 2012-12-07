// float_exponent().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

CL_INLINE sintE CL_INLINE_DECL(float_exponent) (const cl_FF& x)
{
	var uintL uexp = FF_uexp(cl_ffloat_value(x));
	if (uexp==0) { return 0; }
	return (sintL)(uexp - FF_exp_mid);
}

}  // namespace cln
