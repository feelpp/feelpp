// min().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

namespace cln {

const cl_FF min (const cl_FF& x, const cl_FF& y)
{
	return (x <= y ? x : y);
}

}  // namespace cln
