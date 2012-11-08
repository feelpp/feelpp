// round2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_LF_div_t round2 (const cl_LF& x, const cl_LF& y)
{
// Methode:
// (q,r) := round(x/y). Liefere q und x-y*q = y*r.
	var cl_LF_div_t q_r = round2(x/y);
	var cl_I& q = q_r.quotient;
	var cl_LF& r = q_r.remainder;
	return cl_LF_div_t(q,y*r);
}

}  // namespace cln
