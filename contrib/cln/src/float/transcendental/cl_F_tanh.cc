// tanh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

namespace cln {

CL_INLINE const cl_F CL_INLINE_DECL(tanh) (const cl_F& x)
{
// Methode:
// (/ (sinh x) (cosh x))
	var cosh_sinh_t hyp = cosh_sinh(x);
	return The(cl_F)(hyp.sinh) / The(cl_F)(hyp.cosh);
}

}  // namespace cln
