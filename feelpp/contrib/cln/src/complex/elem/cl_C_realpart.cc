// realpart().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"

namespace cln {

const cl_R realpart (const cl_N& x)
{
	if (realp(x)) {
		DeclareType(cl_R,x);
		return x;
	} else {
		DeclareType(cl_C,x);
		return realpart(x);
	}
}

}  // namespace cln
