// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

#include "base/cl_inline.h"
#include "float/ffloat/elem/cl_FF_minusp.cc"
#include "float/ffloat/elem/cl_FF_zerop.cc"

namespace cln {

CL_INLINE2 const cl_FF CL_INLINE2_DECL(signum) (const cl_FF& x)
{
	if (minusp_inline(x)) { return cl_FF_minus1; } // x<0 -> -1.0
	elif (zerop_inline(x)) { return cl_FF_0; } // x=0 -> 0.0
	else { return cl_FF_1; } // x>0 -> +1.0
}

}  // namespace cln
