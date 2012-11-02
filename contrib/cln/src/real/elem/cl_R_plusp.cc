// plusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

namespace cln {

bool plusp (const cl_R& x)
{
	if (minusp(x))
		return false; // x<0 -> nein
	elif (zerop(x))
		return false; // x=0 -> nein
	else
		return true; // sonst ist x>0.
}

}  // namespace cln
