// truncate2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"
#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F_div_t truncate2 (const cl_F& x)
{
	floatcase(x
	,	var cl_SF q = ftruncate(x); return cl_F_div_t(cl_SF_to_I(q),x-q);
	,	var cl_FF q = ftruncate(x); return cl_F_div_t(cl_FF_to_I(q),x-q);
	,	var cl_DF q = ftruncate(x); return cl_F_div_t(cl_DF_to_I(q),x-q);
	,	var cl_LF q = ftruncate(x); return cl_F_div_t(cl_LF_to_I(q),LF_LF_minus_LF(x,q));
	);
}

}  // namespace cln
