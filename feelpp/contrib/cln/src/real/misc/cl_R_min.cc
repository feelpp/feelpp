// min().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

namespace cln {

const cl_R min (const cl_R& x, const cl_R& y)
{
	return (x <= y ? x : y);
}

}  // namespace cln
