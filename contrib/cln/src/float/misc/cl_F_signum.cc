// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

/* Use inline versions of signum(cl_{SF,FF,DF,LF}) functions */
#include "base/cl_inline2.h"
#include "float/sfloat/misc/cl_SF_signum.cc"
#include "float/ffloat/misc/cl_FF_signum.cc"
#include "float/dfloat/misc/cl_DF_signum.cc"
#include "float/lfloat/misc/cl_LF_signum.cc"

namespace cln {

const cl_F CL_FLATTEN signum (const cl_F& x)
{
	floatcase(x
	,	return signum_inline(x);
	,	return signum_inline(x);
	,	return signum_inline(x);
	,	return signum_inline(x);
	);
}

}  // namespace cln
