// float_sign().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline2.h"
#include "float/sfloat/misc/cl_SF_sign.cc"
#include "float/ffloat/misc/cl_FF_sign.cc"
#include "float/dfloat/misc/cl_DF_sign.cc"
#include "float/lfloat/misc/cl_LF_sign.cc"

namespace cln {

const cl_F CL_FLATTEN float_sign (const cl_F& x)
{
// Methode: x>=0 -> Ergebnis 1.0; x<0 -> Ergebnis -1.0
	floatcase(x
	,	return float_sign_inline(x);
	,	return float_sign_inline(x);
	,	return float_sign_inline(x);
	,	return float_sign_inline(x);
	);
}

}  // namespace cln
