// cl_FF_to_DF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"

namespace cln {

const cl_DF cl_FF_to_DF (const cl_FF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	FF_decode(x, { return cl_DF_0; }, sign=,exp=,mant=);
	// Mantisse um 52-23=29 Nullbits erweitern:
	#if (cl_word_size==64)
	return encode_DF(sign,exp,(uint64)mant<<(DF_mant_len-FF_mant_len));
	#else
	return encode_DF(sign,exp,mant>>(32-(DF_mant_len-FF_mant_len)),mant<<(DF_mant_len-FF_mant_len));
	#endif
}

}  // namespace cln
