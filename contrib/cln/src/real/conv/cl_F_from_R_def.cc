// cl_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"

namespace cln {

CL_INLINE const cl_F CL_INLINE_DECL(cl_float) (const cl_R& x)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		return cl_float(x);
	} else {
		DeclareType(cl_F,x);
		return x;
	}
}

}  // namespace cln
