// cl_F_to_float().

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

namespace cln {

float float_approx (const cl_F& x)
{
	floatcase(x
	,	return float_approx(x);
	,	return float_approx(x);
	,	return float_approx(x);
	,	return float_approx(x);
	);
}

}  // namespace cln
