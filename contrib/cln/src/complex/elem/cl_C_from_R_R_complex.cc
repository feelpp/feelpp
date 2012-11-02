// complex_C().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

namespace cln {

const cl_N complex_C (const cl_R& a, const cl_R& b)
{
	return allocate_complex(a,b);
}

}  // namespace cln
