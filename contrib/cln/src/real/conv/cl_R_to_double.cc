// cl_R_to_double().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "float/cl_F.h"
#include "cln/integer.h"
#include "cln/rational.h"
#include "cln/float.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"

#if 0

namespace cln {

double double_approx (const cl_R& x)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		return double_approx(x);
	} else {
		DeclareType(cl_F,x);
		return double_approx(x);
	}
}

}  // namespace cln

#else // fully inlined, faster

#include "rational/cl_RA.h"
#include "integer/cl_I.h"

namespace cln {

double double_approx (const cl_R& x)
{
	realcase6(x
	,	return double_approx(x);
	,	return double_approx(x);
	,	return double_approx(x);
	,	return double_approx(x);
	,	return double_approx(x);
	,	return double_approx(x);
	);
}

}  // namespace cln

#endif
