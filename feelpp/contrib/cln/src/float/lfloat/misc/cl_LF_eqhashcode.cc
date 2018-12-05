// cl_LF equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "base/cl_N.h"
#include "float/lfloat/cl_LF_impl.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

CL_INLINE uint32 CL_INLINE_DECL(equal_hashcode) (const cl_LF& x)
{
	var cl_signean sign;
	var sintL exp;
	var const uintD* MSDptr;
	LF_decode(x, { return 0; }, sign=,exp=,MSDptr=,,);
	#if (intDsize==64)
	var uint32 msd = mspref(MSDptr,0) >> 32;
	#else // (intDsize<=32)
	var uint32 msd = get_32_Dptr(MSDptr);
	#endif
	return equal_hashcode_low(msd,exp,sign);
}

}  // namespace cln
