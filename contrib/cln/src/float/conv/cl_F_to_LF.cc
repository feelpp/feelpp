// cl_F_to_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_LF cl_F_to_LF (const cl_F& x, uintC len)
{
	floatcase(x
	,	return cl_SF_to_LF(x,len);
	,	return cl_FF_to_LF(x,len);
	,	return cl_DF_to_LF(x,len);
	,	return LF_to_LF(x,len);
	);
}

}  // namespace cln
