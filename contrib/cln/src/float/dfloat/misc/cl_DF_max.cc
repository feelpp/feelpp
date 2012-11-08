// max().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

namespace cln {

const cl_DF max (const cl_DF& x, const cl_DF& y)
{
	return (x >= y ? x : y);
}

}  // namespace cln
