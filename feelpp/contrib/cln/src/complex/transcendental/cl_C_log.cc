// log().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"
#include "base/cl_N.h"

namespace cln {

const cl_N log (const cl_N& x)
{
// Methode:
// (complex (log (abs x)) (phase x))
	var cl_R r = abs(x);
	if (zerop(r)) // (abs x) = 0 -> Error
		{ throw division_by_0_exception(); }
	return complex(ln(r),phase(x));
}

}  // namespace cln
