// decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "integer/cl_I.h"

namespace cln {

const decoded_sfloat decode_float (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return decoded_sfloat(SF_0, 0, SF_1); },
		     sign=,exp=,mant=
		 );
	return decoded_sfloat(
		encode_SF(0,0,mant), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_SF(sign,1,bit(SF_mant_len)) // (-1)^s erzeugen
	       );
}

}  // namespace cln
