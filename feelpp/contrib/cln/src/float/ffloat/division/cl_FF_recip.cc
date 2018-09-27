// recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF recip (const cl_FF& x)
{
	return cl_FF_1 / x;
}

}  // namespace cln
