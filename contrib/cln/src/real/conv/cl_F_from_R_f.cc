// cl_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/float.h"
#include "float/cl_F.h"

namespace cln {

const cl_F cl_float (const cl_R& x, float_format_t f)
{
	floatformatcase((uintC)f
	,	return cl_R_to_SF(x);
	,	return cl_R_to_FF(x);
	,	return cl_R_to_DF(x);
	,	return cl_R_to_LF(x,len);
	);
}

}  // namespace cln
