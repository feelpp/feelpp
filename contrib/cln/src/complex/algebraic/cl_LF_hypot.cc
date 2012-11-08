// cl_hypot().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"

/* For inline version of minusp */
#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

ALL_cl_LF_OPERATIONS_SAME_PRECISION()

const cl_LF cl_hypot (const cl_LF& a, const cl_LF& b)
{
//  Zuerst a und b auf gleiche Länge bringen: den längeren runden.
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
 {	Mutable(cl_LF,a);
	Mutable(cl_LF,b);
	{
		var uintC a_len = TheLfloat(a)->len;
		var uintC b_len = TheLfloat(b)->len;
		if (!(a_len == b_len)) {
			if (a_len < b_len)
				b = shorten(b,a_len);
			else
				a = shorten(a,b_len);
		}
	}
	var sintE a_exp;
	var sintE b_exp;
	{
		// Exponenten von a holen:
		var uintE uexp = TheLfloat(a)->expo;
		if (uexp == 0)
			// a=0.0 -> liefere (abs b) :
			return (minusp_inline(b) ? -b : b);
		a_exp = (sintE)(uexp - LF_exp_mid);
	}
	{
		// Exponenten von b holen:
		var uintE uexp = TheLfloat(b)->expo;
		if (uexp == 0)
			// b=0.0 -> liefere (abs a) :
			return (minusp_inline(a) ? -a : a);
		b_exp = (sintE)(uexp - LF_exp_mid);
	}
	// Nun a_exp = float_exponent(a), b_exp = float_exponent(b).
	var sintE e = (a_exp > b_exp ? a_exp : b_exp); // Maximum der Exponenten
	// a und b durch 2^e dividieren:
	var cl_LF na = ((b_exp > a_exp) && ((uintE)(b_exp-a_exp) > (uintE)floor(LF_exp_mid-LF_exp_low-1,2)) ? encode_LF0(TheLfloat(a)->len) : scale_float(a,-e));
	var cl_LF nb = ((a_exp > b_exp) && ((uintE)(a_exp-b_exp) > (uintE)floor(LF_exp_mid-LF_exp_low-1,2)) ? encode_LF0(TheLfloat(b)->len) : scale_float(b,-e));
	// c' := a'*a'+b'*b' berechnen:
	var cl_LF nc = square(na) + square(nb);
	return scale_float(sqrt(nc),e); // c' := sqrt(c'), 2^e*c'
}}

}  // namespace cln
