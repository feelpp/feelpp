// tanh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

namespace cln {

CL_INLINE const cl_R CL_INLINE_DECL(tanh) (const cl_R& x)
{
// Methode:
// (/ (sinh x) (cosh x))
	var cosh_sinh_t hyp = cosh_sinh(x);
	return hyp.sinh / hyp.cosh;
}

}  // namespace cln
