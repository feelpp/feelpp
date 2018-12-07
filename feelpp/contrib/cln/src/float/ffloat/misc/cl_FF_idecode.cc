// integer_decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "float/ffloat/cl_FF.h"
#include "integer/cl_I.h"

#if defined(__mips__) && !defined(__GNUC__)

// Workaround SGI Irix 6.2 C++ 7.0 bug.

#include "float/cl_F.h"
#include "dfloat/cl_DF.h"

namespace cln {

CL_INLINE const cl_idecoded_float CL_INLINE_DECL(integer_decode_float) (const cl_FF& x)
{
	var cl_idecoded_float sem = integer_decode_float(cl_FF_to_DF(x));
	return cl_idecoded_float(sem.mantissa >> (DF_mant_len-FF_mant_len),
	                         sem.exponent + (DF_mant_len-FF_mant_len),
	                         sem.sign
	                        );
}

}  // namespace cln

#else

namespace cln {

CL_INLINE const cl_idecoded_float CL_INLINE_DECL(integer_decode_float) (const cl_FF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	FF_decode(x, { return cl_idecoded_float(0, 0, 1); },
		     sign=,exp=,mant=
		 );
	return cl_idecoded_float(
		(FF_mant_len+1 < cl_value_len
		 ? L_to_FN(mant) // Mantisse als Fixnum
		 : UL_to_I(mant) // oder evtl. als Bignum
		),
		L_to_FN(exp-(FF_mant_len+1)), // e-24 als Fixnum
		(sign>=0 ? cl_I(1) : cl_I(-1)) // (-1)^s erzeugen
	       );
}

}  // namespace cln

#endif
