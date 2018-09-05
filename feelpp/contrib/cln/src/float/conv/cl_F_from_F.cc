// cl_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F cl_float (const cl_F& x, const cl_F& y)
{
	floattypecase(y
	,	return cl_F_to_SF(x);
	,	return cl_F_to_FF(x);
	,	return cl_F_to_DF(x);
	,	return cl_F_to_LF(x,TheLfloat(y)->len);
	);
}

}  // namespace cln
