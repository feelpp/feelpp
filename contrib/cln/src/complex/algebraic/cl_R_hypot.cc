// cl_hypot().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

#include "cln/real.h"
#include "real/cl_R.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"
#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_R cl_hypot (const cl_R& a, const cl_R& b)
{
// Methode:
// Falls a=0: (abs b).
// Falls b=0: (abs a).
// Falls a und b beide rational sind:
//   c:=a*a+b*b, liefere (sqrt c).
// Falls a oder b Floats sind:
//   Falls einer von beiden rational ist, runde ihn zum selben Float-Typ
//     wie der andere und führe das UP durch.
//   Falls beide Floats sind, erweitere auf den genaueren, führe das UP
//     durch und runde wieder auf den ungenaueren.
//   Das Ergebnis ist ein Float >=0.
// UP: [a,b Floats vom selben Typ]
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

	if (rationalp(a)) {
		DeclareType(cl_RA,a);
		if (eq(a,0)) // a=0 -> (abs b)
			return abs(b);
		if (rationalp(b)) {
			DeclareType(cl_RA,b);
			// a,b beide rational
			return sqrt(square(a)+square(b));
		} else {
			DeclareType(cl_F,b);
			// a rational, b Float
			floatcase(b
			,	return cl_hypot(cl_RA_to_SF(a),b);
			,	return cl_hypot(cl_RA_to_FF(a),b);
			,	return cl_hypot(cl_RA_to_DF(a),b);
			,	return cl_hypot(cl_RA_to_LF(a,TheLfloat(b)->len),b);
			);
		}
	} else {
		DeclareType(cl_F,a);
		if (rationalp(b)) {
			DeclareType(cl_RA,b);
			// a Float, b rational
			if (eq(b,0)) // b=0 -> (abs a)
				return abs(a);
			floatcase(a
			,	return cl_hypot(a,cl_RA_to_SF(b));
			,	return cl_hypot(a,cl_RA_to_FF(b));
			,	return cl_hypot(a,cl_RA_to_DF(b));
			,	return cl_hypot(a,cl_RA_to_LF(b,TheLfloat(a)->len));
			);
		} else {
			DeclareType(cl_F,b);
			// a,b Floats
			#ifndef CL_LF_PEDANTIC
			GEN_F_OP2(a,b, cl_hypot, 1, 1, return);
			#else
			GEN_F_OP2(a,b, cl_hypot, 1, 0, return);
			#endif
		}
	}
}

}  // namespace cln
