// ffloor2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "real/division/cl_R_div_t.h"

#if 0 // 2 type dispatches

#include "cln/rational.h"
#include "cln/float.h"
#include "real/division/cl_R_div_t.h"

namespace cln {

const cl_R_fdiv_t ffloor2 (const cl_R& x)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		var cl_RA_div_t q_r = floor2(x);
		var cl_I& q = q_r.quotient;
		var cl_RA& r = q_r.remainder;
		return cl_R_fdiv_t(cl_float(q),r);
	} else {
		DeclareType(cl_F,x);
		return ffloor2(x);
	}
}

}  // namespace cln

#else // 1 type dispatch

#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_R_fdiv_t ffloor2 (const cl_R& x)
{
	realcase6(x
	,	return cl_R_fdiv_t(cl_float(x),0);
	,	var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		var cl_I_div_t q_r = floor2(a,b);
		var cl_I& q = q_r.quotient;
		var cl_I& r = q_r.remainder;
		return cl_R_fdiv_t(cl_float(q),I_I_to_RT(r,b));
	,	var cl_SF q = ffloor(x); return cl_R_fdiv_t(q,x-q);
	,	var cl_FF q = ffloor(x); return cl_R_fdiv_t(q,x-q);
	,	var cl_DF q = ffloor(x); return cl_R_fdiv_t(q,x-q);
	,	var cl_LF q = ffloor(x); return cl_R_fdiv_t(q,LF_LF_minus_LF(x,q));
	);
}

}  // namespace cln

#endif
