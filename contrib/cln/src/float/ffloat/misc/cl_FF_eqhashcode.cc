// cl_FF equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "base/cl_N.h"
#include "float/ffloat/cl_FF.h"

namespace cln {

CL_INLINE uint32 CL_INLINE_DECL(equal_hashcode) (const cl_FF& x)
{
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	FF_decode(x, { return 0; }, sign=,exp=,mant=);
	var uint32 msd = mant << (32-(FF_mant_len+1));
	return equal_hashcode_low(msd,exp,sign);
}

}  // namespace cln
