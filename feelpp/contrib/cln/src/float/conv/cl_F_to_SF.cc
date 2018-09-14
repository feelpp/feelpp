// cl_F_to_SF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

namespace cln {

const cl_SF cl_F_to_SF (const cl_F& x)
{
	floatcase(x
	,	return x;
	,	return cl_FF_to_SF(x);
	,	return cl_DF_to_SF(x);
	,	return cl_LF_to_SF(x);
	);
}

}  // namespace cln
