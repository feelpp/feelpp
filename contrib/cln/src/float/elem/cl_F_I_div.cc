// binary operator /

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

const cl_F operator/ (const cl_F& x, const cl_I& y)
{
	floatcase(x
	, /* SF */	return x / cl_I_to_SF(y);
	, /* FF */	return x / cl_I_to_FF(y);
	, /* DF */	return x / cl_I_to_DF(y);
	, /* LF */	return cl_LF_I_div(x,y);
	);
}

}  // namespace cln
