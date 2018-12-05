// abs().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I abs (const cl_I& x)
{
	// Methode:
	// Bei x<0: (- x), sonst x.
	if (minusp(x))
		return -x;
	else
		return x;
}

}  // namespace cln
