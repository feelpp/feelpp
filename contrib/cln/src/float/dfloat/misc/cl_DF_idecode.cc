// integer_decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"
#include "integer/cl_I.h"

namespace cln {

CL_INLINE const cl_idecoded_float CL_INLINE_DECL(integer_decode_float) (const cl_DF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return cl_idecoded_float(0, 0, 1); },
		     sign=,exp=,mant=
		 );
	return cl_idecoded_float(
		Q_to_I(mant), // Mantisse (>0, <2^53) als Bignum
		L_to_FN(exp-(DF_mant_len+1)), // e-53 als Fixnum
		(sign>=0 ? cl_I(1) : cl_I(-1)) // (-1)^s erzeugen
	       );
#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return cl_idecoded_float(0, 0, 1); },
		      sign=,exp=,manthi=,mantlo=
		  );
	return cl_idecoded_float(
		L2_to_I(manthi,mantlo), // Mantisse (>0, <2^53) als Bignum
		L_to_FN(exp-(DF_mant_len+1)), // e als Fixnum
		(sign>=0 ? cl_I(1) : cl_I(-1)) // (-1)^s erzeugen
	       );
#endif
}

}  // namespace cln
