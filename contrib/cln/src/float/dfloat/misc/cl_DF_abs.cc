// abs().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "base/cl_inline.h"
#include "float/dfloat/elem/cl_DF_minusp.cc"

namespace cln {

const cl_DF CL_FLATTEN abs (const cl_DF& x)
{
// x<0 -> (- x), x>=0 -> x
	if (minusp(x)) return -x; else return x;
}

}  // namespace cln
