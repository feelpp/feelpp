// decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"
#include "integer/cl_I.h"

namespace cln {

const decoded_lfloat decode_float (const cl_LF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintE exp;
	var uintC mantlen;
	var const uintD* mantMSDptr;
	LF_decode(x, { return decoded_lfloat(x, 0, encode_LF1(mantlen)); },
		     sign=,exp=,mantMSDptr=,mantlen=,);
	return decoded_lfloat(
		encode_LFu(0,0+LF_exp_mid,mantMSDptr,mantlen), // (-1)^0 * 2^0 * m erzeugen
		E_to_I(exp), // e als Fixnum
		encode_LF1s(sign,mantlen) // (-1)^s erzeugen
	       );
}

}  // namespace cln
