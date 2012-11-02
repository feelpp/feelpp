// fceiling().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"

#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

CL_INLINE2 const cl_LF CL_INLINE2_DECL(fceiling) (const cl_LF& x)
{
	if (minusp_inline(x))
		return ftruncate(x);
	else
		return futruncate(x);
}

}  // namespace cln
