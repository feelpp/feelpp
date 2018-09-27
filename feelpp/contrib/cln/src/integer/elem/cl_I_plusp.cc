// plusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#define minusp inline_minusp
#define zerop inline_zerop
#include "integer/cl_I.h"
#undef zerop
#undef minusp

namespace cln {

bool plusp (const cl_I& x)
{
	if (inline_minusp(x))
		return false; // x<0 -> nein
	elif (inline_zerop(x))
		return false; // x=0 -> nein
	else
		return true; // sonst ist x>0.
}

}  // namespace cln
