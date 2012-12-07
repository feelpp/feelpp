// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"
#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_R cl_F_RA_mul (const cl_F& x, const cl_RA& y)
{
	if (eq(y,0)) { return 0; } // x * 0 = exakte 0
	floatcase(x
	, /* SF */	return x * cl_RA_to_SF(y);
	, /* FF */	return x * cl_RA_to_FF(y);
	, /* DF */	return x * cl_RA_to_DF(y);
	, /* LF */	return cl_LF_RA_mul(x,y);
	);
}

}  // namespace cln
