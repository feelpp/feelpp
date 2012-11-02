// cl_SF_to_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"

namespace cln {

float float_approx (const cl_SF& x)
{
	// x entpacken:
	var cl_signean sign;
	var sintL exp;
	var uint32 mant;
	SF_decode(x, { return 0.0; }, sign=,exp=,mant=);
	// Mantisse um 23-16=7 Bits nach links schieben:
	union { ffloat eksplicit; float machine_float; } u;
	if (((sintL)(SF_exp_high-SF_exp_mid) > (sintL)(FF_exp_high-FF_exp_mid))
	    && (exp > (sintL)(FF_exp_high-FF_exp_mid)))
	  { u.eksplicit = make_FF_word(sign,bit(FF_exp_len)-1,0); } // Infinity
	else
	if (((sintL)(SF_exp_low-SF_exp_mid) < (sintL)(FF_exp_low-FF_exp_mid))
	    && (exp < (sintL)(FF_exp_low-FF_exp_mid)))
	  { u.eksplicit = make_FF_word(sign,0,0); } // 0.0
	else
	  { u.eksplicit = make_FF_word(sign,exp+FF_exp_mid,mant<<(FF_mant_len-SF_mant_len)); }
	return u.machine_float;
}

}  // namespace cln
