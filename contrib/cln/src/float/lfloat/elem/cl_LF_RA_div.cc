// cl_LF_RA_div().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/lfloat.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"

namespace cln {

const cl_LF cl_LF_RA_div (const cl_LF& x, const cl_RA& y)
{
// Method:
// Write y = u/v. Return (x*v)/u.
	if (integerp(y)) {
		DeclareType(cl_I,y);
		return cl_LF_I_div(x,y);
	} else {
		DeclareType(cl_RT,y);
		var const cl_I& u = TheRatio(y)->numerator; // u /= 0
		var const cl_I& v = TheRatio(y)->denominator; // v /= 0
		return cl_LF_I_div(The(cl_LF)(cl_LF_I_mul(x,v)),u);
	}
}

}  // namespace cln
