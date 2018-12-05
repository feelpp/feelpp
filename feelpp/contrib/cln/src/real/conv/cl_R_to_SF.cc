// cl_R_to_SF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"

#if 0

namespace cln {

const cl_SF cl_R_to_SF (const cl_R& x)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		return cl_RA_to_SF(x);
	} else {
		DeclareType(cl_F,x);
		return cl_F_to_SF(x);
	}
}

}  // namespace cln

#else // fully inlined, faster

#include "rational/cl_RA.h"
#include "integer/cl_I.h"

namespace cln {

const cl_SF cl_R_to_SF (const cl_R& x)
{
	realcase6(x
	,	return cl_I_to_SF(x);
	,	return cl_RA_to_SF(x);
	,	return x;
	,	return cl_FF_to_SF(x);
	,	return cl_DF_to_SF(x);
	,	return cl_LF_to_SF(x);
	);
}

}  // namespace cln

#endif
