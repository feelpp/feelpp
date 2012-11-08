// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

#include "base/cl_inline.h"
#include "float/dfloat/elem/cl_DF_minusp.cc"
#include "float/dfloat/elem/cl_DF_zerop.cc"

namespace cln {

CL_INLINE2 const cl_DF CL_INLINE2_DECL(signum) (const cl_DF& x)
{
	if (minusp_inline(x)) { return cl_DF_minus1; } // x<0 -> -1.0
	elif (zerop_inline(x)) { return cl_DF_0; } // x=0 -> 0.0
	else { return cl_DF_1; } // x>0 -> +1.0
}

}  // namespace cln
