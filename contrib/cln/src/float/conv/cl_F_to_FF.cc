// cl_F_to_FF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

namespace cln {

const cl_FF cl_F_to_FF (const cl_F& x)
{
	floatcase(x
	,	return cl_SF_to_FF(x);
	,	return x;
	,	return cl_DF_to_FF(x);
	,	return cl_LF_to_FF(x);
	);
}

}  // namespace cln
