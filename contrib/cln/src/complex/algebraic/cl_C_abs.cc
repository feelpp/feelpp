// abs().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

#include "base/cl_inline.h"
#include "complex/algebraic/cl_C_abs_aux.cc"

namespace cln {

const cl_R abs (const cl_N& x)
{
// Methode:
// Falls x reell: klar
// Falls x=a+bi: sqrt(a^2+b^2)
	if (realp(x)) {
		DeclareType(cl_R,x);
		return abs(x);
	} else {
		DeclareType(cl_C,x);
		return abs_inline(x);
	}
}

}  // namespace cln
