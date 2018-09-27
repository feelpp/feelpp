// recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

namespace cln {

const cl_SF recip (const cl_SF& x)
{
	return SF_1 / x;
}

}  // namespace cln
