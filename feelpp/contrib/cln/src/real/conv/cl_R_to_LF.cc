// cl_F_to_LF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "float/cl_F.h"
#include "float/lfloat/cl_LF.h"

#if 0

namespace cln {

const cl_LF cl_R_to_LF (const cl_R& x, uintC len)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		return cl_RA_to_LF(x,len);
	} else {
		DeclareType(cl_F,x);
		return cl_F_to_LF(x,len);
	}
}

}  // namespace cln

#else // fully inlined, faster

#include "rational/cl_RA.h"
#include "integer/cl_I.h"

namespace cln {

const cl_LF cl_R_to_LF (const cl_R& x, uintC len)
{
	realcase6(x
	,	return cl_I_to_LF(x,len);
	,	return cl_RA_to_LF(x,len);
	,	return cl_SF_to_LF(x,len);
	,	return cl_FF_to_LF(x,len);
	,	return cl_DF_to_LF(x,len);
	,	return LF_to_LF(x,len);
	);
}

}  // namespace cln

#endif
