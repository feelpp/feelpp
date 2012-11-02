// cl_hypot().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

#include "cln/ffloat.h"
#include "float/ffloat/cl_FF.h"

/* For inline version of minusp */
#include "base/cl_inline.h"
#include "float/ffloat/elem/cl_FF_minusp.cc"

namespace cln {

const cl_FF cl_hypot (const cl_FF& a, const cl_FF& b)
{
//  a=0.0 -> liefere abs(b).
//  b=0.0 -> liefere abs(a).
//  e:=max(exponent(a),exponent(b)).
//  a':=a/2^e bzw. 0.0 bei Underflowmöglichkeit (beim Skalieren a':=a/2^e
//      oder beim Quadrieren a'*a':  2*(e-exponent(a))>exp_mid-exp_low-1
//      d.h. exponent(b)-exponent(a)>floor((exp_mid-exp_low-1)/2) ).
//  b':=b/2^e bzw. 0.0 bei Underflowmöglichkeit (beim Skalieren b':=b/2^e
//      oder beim Quadrieren b'*b':  2*(e-exponent(b))>exp_mid-exp_low-1
//      d.h. exponent(a)-exponent(b)>floor((exp_mid-exp_low-1)/2) ).
//  c':=a'*a'+b'*b', c':=sqrt(c'), liefere 2^e*c'.
	var sintL a_exp;
	var sintL b_exp;
	{
		// Exponenten von a holen:
		var uintL uexp = FF_uexp(cl_ffloat_value(a));
		if (uexp == 0)
			// a=0.0 -> liefere (abs b) :
			return (minusp_inline(b) ? -b : b);
		a_exp = (sintL)(uexp - FF_exp_mid);
	}
	{
		// Exponenten von b holen:
		var uintL uexp = FF_uexp(cl_ffloat_value(b));
		if (uexp == 0)
			// b=0.0 -> liefere (abs a) :
			return (minusp_inline(a) ? -a : a);
		b_exp = (sintL)(uexp - FF_exp_mid);
	}
	// Nun a_exp = float_exponent(a), b_exp = float_exponent(b).
	var sintL e = (a_exp > b_exp ? a_exp : b_exp); // Maximum der Exponenten
	// a und b durch 2^e dividieren:
	var cl_FF na = (b_exp-a_exp > floor(FF_exp_mid-FF_exp_low-1,2) ? cl_FF_0 : scale_float(a,-e));
	var cl_FF nb = (a_exp-b_exp > floor(FF_exp_mid-FF_exp_low-1,2) ? cl_FF_0 : scale_float(b,-e));
	// c' := a'*a'+b'*b' berechnen:
	var cl_FF nc = square(na) + square(nb);
	return scale_float(sqrt(nc),e); // c' := sqrt(c'), 2^e*c'
}

}  // namespace cln
