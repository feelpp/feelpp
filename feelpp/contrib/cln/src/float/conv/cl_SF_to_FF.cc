// cl_SF_to_FF().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"

namespace cln {

const cl_FF cl_SF_to_FF (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return cl_FF_0; }, sign=,exp=,mant=);
	// Mantisse um 23-16=7 Bits nach links schieben:
	return encode_FF(sign,exp,mant<<(FF_mant_len-SF_mant_len));
}

}  // namespace cln
