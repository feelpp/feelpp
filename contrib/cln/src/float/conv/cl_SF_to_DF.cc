// cl_SF_to_DF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/dfloat/cl_DF.h"

namespace cln {

const cl_DF cl_SF_to_DF (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return cl_DF_0; }, sign=,exp=,mant=);
	// Mantisse um 52-16=36 Nullbits erweitern:
	#if (cl_word_size==64)
	return encode_DF(sign,exp,(uint64)mant<<(DF_mant_len-SF_mant_len));
	#else
	return encode_DF(sign,exp,mant<<(DF_mant_len-SF_mant_len-32),0);
	#endif
}

}  // namespace cln
