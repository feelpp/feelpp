// min().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

namespace cln {

const cl_F min (const cl_F& x, const cl_F& y)
{
	return (x <= y ? x : y);
}

}  // namespace cln
