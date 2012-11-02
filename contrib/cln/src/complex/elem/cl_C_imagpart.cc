// imagpart().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"

namespace cln {

const cl_R imagpart (const cl_N& x)
{
	if (realp(x))
		return 0;
	else {
		DeclareType(cl_C,x);
		return imagpart(x);
	}
}

}  // namespace cln
