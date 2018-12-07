// cl_DF equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "base/cl_N.h"
#include "float/dfloat/cl_DF.h"

namespace cln {

CL_INLINE uint32 CL_INLINE_DECL(equal_hashcode) (const cl_DF& x)
{
	var cl_signean sign;
	var sintL exp;
#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return 0; }, sign=,exp=,mant=);
	var uint32 msd = mant >> ((DF_mant_len+1)-32);
#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return 0; }, sign=,exp=,manthi=,mantlo=);
	var uint32 msd = (manthi << (64-(DF_mant_len+1)))
			 | (mantlo >> ((DF_mant_len+1)-32));
#endif
	return equal_hashcode_low(msd,exp,sign);
}

}  // namespace cln
