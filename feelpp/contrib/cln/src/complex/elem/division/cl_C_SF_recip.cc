// cl_C_recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

#include "cln/sfloat.h"
#include "float/sfloat/cl_SF.h"

namespace cln {

const cl_C_SF cl_C_recip (const cl_SF& a, const cl_SF& b)
{
//  a=0.0 -> liefere die Komponenten a=0.0 und -1/b.
//  b=0.0 -> liefere die Komponenten 1/a und b=0.0.
//  e:=max(exponent(a),exponent(b)).
//  a':=a/2^e bzw. 0.0 bei Underflowmöglichkeit (beim Skalieren a':=a/2^e
//      oder beim Quadrieren a'*a':  2*(e-exponent(a))>exp_mid-exp_low-1
//      d.h. exponent(b)-exponent(a)>floor((exp_mid-exp_low-1)/2) ).
//  b':=b/2^e bzw. 0.0 bei Underflowmöglichkeit (beim Skalieren b':=b/2^e
//      oder beim Quadrieren b'*b':  2*(e-exponent(b))>exp_mid-exp_low-1
//      d.h. exponent(a)-exponent(b)>floor((exp_mid-exp_low-1)/2) ).
//  c':=a'*a'+b'*b',
//  liefere die beiden Komponenten 2^(-e)*a'/c' und -2^(-e)*b'/c'.
	var sintL a_exp;
	var sintL b_exp;
	{
		// Exponenten von a holen:
		var uintL uexp = SF_uexp(a);
		if (uexp == 0)
			// a=0.0 -> liefere (complex a (- (/ b))) :
			return cl_C_SF(a,-recip(b));
		a_exp = (sintL)(uexp - SF_exp_mid);
	}
	{
		// Exponenten von b holen:
		var uintL uexp = SF_uexp(b);
		if (uexp == 0)
			// b=0.0 -> liefere (complex (/ a) b) :
			return cl_C_SF(recip(a),b);
		b_exp = (sintL)(uexp - SF_exp_mid);
	}
	// Nun a_exp = float_exponent(a), b_exp = float_exponent(b).
	var sintL e = (a_exp > b_exp ? a_exp : b_exp); // Maximum der Exponenten
	// a und b durch 2^e dividieren:
	var cl_SF na = (b_exp-a_exp > floor(SF_exp_mid-SF_exp_low-1,2) ? SF_0 : scale_float(a,-e));
	var cl_SF nb = (a_exp-b_exp > floor(SF_exp_mid-SF_exp_low-1,2) ? SF_0 : scale_float(b,-e));
	// c' := a'*a'+b'*b' berechnen:
	var cl_SF nc = square(na) + square(nb);
	// 2^(-e)*a'/c' + i * -2^(-e)*b'/c'
	return cl_C_SF(scale_float(na/nc,-e), scale_float(-(nb/nc),-e));
}

}  // namespace cln
