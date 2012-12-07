// decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "integer/cl_I.h"

namespace cln {

const decoded_dfloat decode_float (const cl_DF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return decoded_dfloat(cl_DF_0, 0, cl_DF_1); },
		     sign=,exp=,mant=
		 );
	return decoded_dfloat(
		encode_DF(0,0,mant), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_DF(sign,1,bit(DF_mant_len)) // (-1)^s erzeugen
	       );
#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return decoded_dfloat(cl_DF_0, 0, cl_DF_1); },
		      sign=,exp=,manthi=,mantlo=
		  );
	return decoded_dfloat(
		encode_DF(0,0,manthi,mantlo), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_DF(sign,1,bit(DF_mant_len-32),0) // (-1)^s erzeugen
	       );
#endif
}

}  // namespace cln
