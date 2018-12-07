// minusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#define minusp inline_minusp
#include "integer/cl_I.h"
#undef minusp

namespace cln {

bool minusp (const cl_I& x)
{
	return inline_minusp(x);
}

}  // namespace cln
