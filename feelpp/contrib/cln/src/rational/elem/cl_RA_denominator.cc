// denominator().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"

namespace cln {

const cl_I denominator (const cl_RA& r)
{
	if (integerp(r))
		return 1;
	else
		return TheRatio(r)->denominator;
}

}  // namespace cln
