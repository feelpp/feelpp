// minusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

CL_INLINE bool CL_INLINE_DECL(minusp) (const cl_FF& x)
{
	return (sint32)cl_ffloat_value(x) < 0;
}

}  // namespace cln
