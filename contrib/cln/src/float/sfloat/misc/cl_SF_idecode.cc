// integer_decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "integer/cl_I.h"

namespace cln {

CL_INLINE const cl_idecoded_float CL_INLINE_DECL(integer_decode_float) (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return cl_idecoded_float(0, 0, 1); },
		     sign=,exp=,mant=
		 );
	return cl_idecoded_float(
		L_to_FN(mant), // Mantisse als Fixnum (>0, <2^17)
		L_to_FN(exp-(SF_mant_len+1)), // e-17 als Fixnum
		(sign>=0 ? cl_I(1) : cl_I(-1)) // (-1)^s erzeugen
	       );
}

}  // namespace cln
