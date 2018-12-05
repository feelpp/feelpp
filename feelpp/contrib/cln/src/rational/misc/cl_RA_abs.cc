// abs().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA abs (const cl_RA& r)
{
	// Methode:
	// Bei r<0: (- r), sonst r.
	if (minusp(r))
		return -r;
	else
		return r;
}

}  // namespace cln
