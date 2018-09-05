// minusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#define minusp minusp_inline
#include "rational/cl_RA.h"
#undef minusp

namespace cln {

bool CL_FLATTEN minusp (const cl_RA& x)
{
	return minusp_inline(x);
}

}  // namespace cln
