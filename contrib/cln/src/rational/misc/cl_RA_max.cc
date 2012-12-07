// max().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

namespace cln {

const cl_RA max (const cl_RA& x, const cl_RA& y)
{
	return (x >= y ? x : y);
}

}  // namespace cln
