// recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
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

ALL_cl_LF_OPERATIONS_SAME_PRECISION()

// for GEN_F_OP2:
#define NOMAP2(F,EXPR)  \
  cl_C_##F _tmp = EXPR;							\
  return complex_C(_tmp.realpart,_tmp.imagpart);
#define MAP2(F,FN,EXPR)  \
  cl_C_##F _tmp = EXPR;							\
  return complex_C(FN(_tmp.realpart),FN(_tmp.imagpart))

const cl_N recip (const cl_N& x)
{
// Methode:
// Falls x reell, klar.
// Falls x=a+bi:
//    Falls a=0: 0+(-1/b)i.
//    Falls a und b beide rational sind:
//      c:=a*a+b*b, c:=1/c, liefere a*c+(-b*c)i.
//    Falls a oder b Floats sind:
//      Falls einer von beiden rational ist, runde ihn zum selben Float-Typ
//        wie der andere und führe das UP durch.
//      Falls beide Floats sind, erweitere auf den genaueren, führe das UP
//        durch und runde wieder auf den ungenaueren.
//      Das Ergebnis ist eine komplexe Zahl, da beide Komponenten Floats sind.
// UP: [a,b Floats vom selben Typ]
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

	if (realp(x)) {
		DeclareType(cl_R,x);
		return recip(x);
	} else
    {
	DeclareType(cl_C,x);
	var const cl_R& a = realpart(x);
	var const cl_R& b = imagpart(x);
	// x = a+bi
	if (rationalp(a)) {
		DeclareType(cl_RA,a);
		if (eq(a,0))
			// (complex 0 (- (/ b)))
			return complex_C(0,-recip(b));
		if (rationalp(b)) {
			DeclareType(cl_RA,b);
			// a,b beide rational
			var cl_RA c = recip(square(a)+square(b));
			return complex_C(a*c,-b*c);
		} else {
			DeclareType(cl_F,b);
			// a rational, b Float
			floatcase(b
			,	return complex_C(cl_C_recip(cl_RA_to_SF(a),b));
			,	return complex_C(cl_C_recip(cl_RA_to_FF(a),b));
			,	return complex_C(cl_C_recip(cl_RA_to_DF(a),b));
			,	return complex_C(cl_C_recip(cl_RA_to_LF(a,TheLfloat(b)->len),b));
			);
		}
	} else {
		DeclareType(cl_F,a);
		if (rationalp(b)) {
			DeclareType(cl_RA,b);
			// a Float, b rational
			floatcase(a
			,	return complex_C(cl_C_recip(a,cl_RA_to_SF(b)));
			,	return complex_C(cl_C_recip(a,cl_RA_to_FF(b)));
			,	return complex_C(cl_C_recip(a,cl_RA_to_DF(b)));
			,	return complex_C(cl_C_recip(a,cl_RA_to_LF(b,TheLfloat(a)->len)));
			);
		} else {
			DeclareType(cl_F,b);
			// a,b Floats
			#ifndef CL_LF_PEDANTIC
			GEN_F_OP2(a,b,cl_C_recip,2,1,); // uses NOMAP2, MAP2.
			#else
			GEN_F_OP2(a,b,cl_C_recip,2,0,); // uses NOMAP2, MAP2.
			#endif
		}
	}
    }
}

}  // namespace cln
