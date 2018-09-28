// zerop().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#define zerop inline_zerop
#include "rational/cl_RA.h"
#undef zerop

namespace cln {

bool zerop (const cl_RA& x)
{
	return inline_zerop(x);
}

}  // namespace cln
