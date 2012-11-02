// abs().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_minusp.cc"

namespace cln {

const cl_SF CL_FLATTEN abs (const cl_SF& x)
{
// x<0 -> (- x), x>=0 -> x
	if (minusp_inline(x)) return -x; else return x;
}

}  // namespace cln
