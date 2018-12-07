// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "integer/cl_I.h"


namespace cln {

CL_INLINE const cl_RA CL_INLINE_DECL(signum) (const cl_RA& x)
{
	if (minusp(x)) { return -1; } // x<0 -> -1
	elif (zerop(x)) { return 0; } // x=0 -> 0
	else { return 1; } // x>0 -> +1
}

}  // namespace cln
