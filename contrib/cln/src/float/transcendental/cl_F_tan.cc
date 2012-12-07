// tan().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

namespace cln {

CL_INLINE const cl_F CL_INLINE_DECL(tan) (const cl_F& x)
{
// Methode:
// (/ (sin x) (cos x))
	var cos_sin_t trig = cos_sin(x);
	return The(cl_F)(trig.sin) / The(cl_F)(trig.cos);
}

}  // namespace cln
