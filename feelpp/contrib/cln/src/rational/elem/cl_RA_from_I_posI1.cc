// I_I_to_RT().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "rational/cl_RA.h"


// Implementation.

namespace cln {

const cl_RA I_I_to_RT (const cl_I& a, const cl_I& b)
{
	return allocate_ratio(a,b);
}

}  // namespace cln
