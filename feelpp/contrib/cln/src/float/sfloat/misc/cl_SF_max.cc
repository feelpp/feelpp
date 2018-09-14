// max().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

namespace cln {

const cl_SF max (const cl_SF& x, const cl_SF& y)
{
	return (x >= y ? x : y);
}

}  // namespace cln
