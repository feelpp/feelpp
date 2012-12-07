// cl_SF equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "base/cl_N.h"
#include "float/sfloat/cl_SF.h"

namespace cln {

CL_INLINE uint32 CL_INLINE_DECL(equal_hashcode) (const cl_SF& x)
{
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return 0; }, sign=,exp=,mant=);
	var uint32 msd = mant << (32-(SF_mant_len+1));
	return equal_hashcode_low(msd,exp,sign);
}

}  // namespace cln
