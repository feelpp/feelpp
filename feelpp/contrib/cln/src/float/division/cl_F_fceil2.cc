// fceiling2().

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
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F_fdiv_t fceiling2 (const cl_F& x)
{
#if 0 // 3 type dispatches
	var cl_F q = fceiling(x);
	return cl_F_fdiv_t(q,x-q);
#else // 1 type dispatch
	floatcase(x
	,	var cl_SF q = fceiling(x); return cl_F_fdiv_t(q,x-q);
	,	var cl_FF q = fceiling(x); return cl_F_fdiv_t(q,x-q);
	,	var cl_DF q = fceiling(x); return cl_F_fdiv_t(q,x-q);
	,	var cl_LF q = fceiling(x); return cl_F_fdiv_t(q,LF_LF_minus_LF(x,q));
	);
#endif
}

}  // namespace cln
