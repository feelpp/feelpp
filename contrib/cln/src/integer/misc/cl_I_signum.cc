// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

CL_INLINE const cl_I CL_INLINE_DECL(signum) (const cl_I& x)
{
	if (minusp(x)) { return -1; } // x<0 -> -1
	elif (zerop(x)) { return 0; } // x=0 -> 0
	else { return 1; } // x>0 -> +1
}

}  // namespace cln
