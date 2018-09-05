// tan().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

namespace cln {

CL_INLINE const cl_R CL_INLINE_DECL(tan) (const cl_R& x)
{
// Methode:
// (/ (sin x) (cos x))
	var cos_sin_t trig = cos_sin(x);
	return trig.sin / trig.cos;
}

}  // namespace cln
