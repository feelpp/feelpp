// conjugate().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N conjugate (const cl_N& x)
{
	if (realp(x))
		return x;
	else {
		DeclareType(cl_C,x);
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		// Vorzeichenwechsel beim Imaginärteil
		return complex_C(a,-b);
	}
}

}  // namespace cln
