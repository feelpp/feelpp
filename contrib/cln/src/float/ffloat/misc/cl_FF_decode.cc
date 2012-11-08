// decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "integer/cl_I.h"

namespace cln {

const decoded_ffloat decode_float (const cl_FF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	FF_decode(x, { return decoded_ffloat(cl_FF_0, 0, cl_FF_1); },
		     sign=,exp=,mant=
		 );
	return decoded_ffloat(
		encode_FF(0,0,mant), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_FF(sign,1,bit(FF_mant_len)) // (-1)^s erzeugen
	       );
}

}  // namespace cln
