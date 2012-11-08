// fround().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"

namespace cln {

const cl_F fround (const cl_R& x)
{
	realcase6(x
	,	return cl_float(x);
	,	var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		return cl_float(round1(a,b));
	,	return fround(x);
	,	return fround(x);
	,	return fround(x);
	,	return fround(x);
	);
}

}  // namespace cln
