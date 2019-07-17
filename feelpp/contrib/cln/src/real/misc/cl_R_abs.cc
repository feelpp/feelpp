// abs().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

namespace cln {

const cl_R abs (const cl_R& x)
{
	// Methode:
	// Bei x<0: (- x), sonst x.
	if (minusp(x))
		return -x;
	else
		return x;
}

}  // namespace cln
