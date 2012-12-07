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

const cl_R operator/ (const cl_RA& x, const cl_F& y)
{
	if (eq(x,0)) { return 0; }
	floatcase(y
	, /* SF */	if (integerp(x)) {
				DeclareType(cl_I,x);
				return cl_I_to_SF(x) / y;
			} else
				return cl_RA_to_SF(x) / y;
	, /* FF */	if (integerp(x)) {
				DeclareType(cl_I,x);
				return cl_I_to_FF(x) / y;
			} else
				return cl_RA_to_FF(x) / y;
	, /* DF */	if (integerp(x)) {
				DeclareType(cl_I,x);
				return cl_I_to_DF(x) / y;
			} else
				return cl_RA_to_DF(x) / y;
	, /* LF */	if (integerp(x)) {
				DeclareType(cl_I,x);
				return cl_I_LF_div(x,y);
			} else
				return cl_RA_LF_div(x,y);
	);
}

}  // namespace cln
