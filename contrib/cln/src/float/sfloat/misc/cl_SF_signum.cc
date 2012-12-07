// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

/* Use inline versions of minusp and zerop */
#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_minusp.cc"
#include "float/sfloat/elem/cl_SF_zerop.cc"

namespace cln {

CL_INLINE2 const cl_SF CL_INLINE2_DECL(signum) (const cl_SF& x)
{
	if (minusp_inline(x)) { return SF_minus1; } // x<0 -> -1.0
	elif (zerop_inline(x)) { return SF_0; } // x=0 -> 0.0
	else { return SF_1; } // x>0 -> +1.0
}

}  // namespace cln
