// fceiling().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

#include "base/cl_inline.h"
#include "float/ffloat/elem/cl_FF_minusp.cc"

namespace cln {

const cl_FF CL_FLATTEN fceiling (const cl_FF& x)
{
	if (minusp_inline(x))
		return ftruncate(x);
	else
		return futruncate(x);
}

}  // namespace cln
