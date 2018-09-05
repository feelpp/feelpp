// truncate2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

const cl_DF_div_t truncate2 (const cl_DF& x, const cl_DF& y)
{
// Methode:
// (q,r) := truncate(x/y). Liefere q und x-y*q = y*r.
	var cl_DF_div_t q_r = truncate2(x/y);
	var cl_I& q = q_r.quotient;
	var cl_DF& r = q_r.remainder;
	return cl_DF_div_t(q,y*r);
}

}  // namespace cln
