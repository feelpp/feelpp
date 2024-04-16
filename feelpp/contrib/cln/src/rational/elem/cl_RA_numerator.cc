// numerator().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"

namespace cln {

const cl_I numerator (const cl_RA& r)
{
	if (integerp(r)) {
		DeclareType(cl_I,r);
		return r;
	} else
		return TheRatio(r)->numerator;
}

}  // namespace cln
