// plusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#define minusp inline_minusp
#define zerop inline_zerop
#include "rational/cl_RA.h"
#undef zerop
#undef minusp

namespace cln {

bool plusp (const cl_RA& x)
{
	if (inline_minusp(x))
		return false; // x<0 -> nein
	elif (inline_zerop(x))
		return false; // x=0 -> nein
	else
		return true; // sonst ist x>0.
}

}  // namespace cln
