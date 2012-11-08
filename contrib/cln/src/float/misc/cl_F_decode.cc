// decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"

namespace cln {

inline const decoded_float decode_float (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return decoded_float(SF_0, 0, SF_1); },
		     sign=,exp=,mant=
		 );
	return decoded_float(
		encode_SF(0,0,mant), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_SF(sign,1,bit(SF_mant_len)) // (-1)^s erzeugen
	       );
}

inline const decoded_float decode_float (const cl_FF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	FF_decode(x, { return decoded_float(cl_FF_0, 0, cl_FF_1); },
		     sign=,exp=,mant=
		 );
	return decoded_float(
		encode_FF(0,0,mant), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_FF(sign,1,bit(FF_mant_len)) // (-1)^s erzeugen
	       );
}

inline const decoded_float decode_float (const cl_DF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
#if (cl_word_size==64)
	var uint64 mant;
	DF_decode(x, { return decoded_float(cl_DF_0, 0, cl_DF_1); },
		     sign=,exp=,mant=
		 );
	return decoded_float(
		encode_DF(0,0,mant), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_DF(sign,1,bit(DF_mant_len)) // (-1)^s erzeugen
	       );
#else
	var uint32 manthi;
	var uint32 mantlo;
	DF_decode2(x, { return decoded_float(cl_DF_0, 0, cl_DF_1); },
		      sign=,exp=,manthi=,mantlo=
		  );
	return decoded_float(
		encode_DF(0,0,manthi,mantlo), // (-1)^0 * 2^0 * m erzeugen
		L_to_FN(exp), // e als Fixnum
		encode_DF(sign,1,bit(DF_mant_len-32),0) // (-1)^s erzeugen
	       );
#endif
}

inline const decoded_float decode_float (const cl_LF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintE exp;
	var uintC mantlen;
	var const uintD* mantMSDptr;
	LF_decode(x, { return decoded_float(x, 0, encode_LF1(mantlen)); },
		     sign=,exp=,mantMSDptr=,mantlen=,);
	return decoded_float(
		encode_LFu(0,0+LF_exp_mid,mantMSDptr,mantlen), // (-1)^0 * 2^0 * m erzeugen
		E_to_I(exp), // e als Fixnum
		encode_LF1s(sign,mantlen) // (-1)^s erzeugen
	       );
}

const decoded_float decode_float (const cl_F& x)
{
	floatcase(x
	,	return decode_float(x);
	,	return decode_float(x);
	,	return decode_float(x);
	,	return decode_float(x);
	);
}

}  // namespace cln
