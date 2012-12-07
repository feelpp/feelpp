// floor2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF_div_t floor2 (const cl_FF& x, const cl_FF& y)
{
// Methode:
// (q,r) := floor(x/y). Liefere q und x-y*q = y*r.
	var cl_FF_div_t q_r = floor2(x/y);
	var cl_I& q = q_r.quotient;
	var cl_FF& r = q_r.remainder;
	return cl_FF_div_t(q,y*r);
}

}  // namespace cln
