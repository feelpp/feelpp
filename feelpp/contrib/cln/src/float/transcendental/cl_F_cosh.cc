// cosh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/transcendental/cl_F_tran.h"
#include "float/cl_F.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F cosh (const cl_F& x)
{
// Methode:
// Genauigkeit erhöhen,
// e := Exponent aus (decode-float x), d := (float-digits x)
// falls x=0.0 oder e<=(1-d)/2 liefere 1.0
//   (denn bei e<=(1-d)/2 ist 1 <= cosh(x) = 1+x^2/2+... < 1+2^(-d),
//    also ist cosh(x), auf d Bits gerundet, gleich 1.0).
// falls e<0:
//   y := x/2 = (scale-float x -1), (sinh(y)/y)^2 errechnen,
//   cosh(x) = 1+x*y*(sinh(y)/y)^2 errechnen.
// falls e>=0: y:=exp(x) errechnen, (scale-float (+ y (/ y)) -1) bilden.

	var sintE e = float_exponent(x);
	if (e < 0) { // Exponent e abtesten
		// e<0
		if (zerop(x))
			return cl_float(1,x);
		var uintC d = float_digits(x);
		if (e <= (1-(sintC)d)>>1) // e <= (1-d)/2 <==> e <= -ceiling((d-1)/2) ?
			return cl_float(1,x); // ja -> 1.0 als Ergebnis
		// Rechengenauigkeit erhöhen
		if (longfloatp(x)) {
			DeclareType(cl_LF,x);
			#if 0
			if (TheLfloat(x)->len >= infty) {
				var cl_LF xx = extend(x,TheLfloat(x)->len+1);
				var cl_LF_cosh_sinh_t hyp = cl_coshsinh_ratseries(xx);
				return cl_float(hyp.cosh,x);
			} else
			#endif
			if (TheLfloat(x)->len >= 600) {
				// verwende exp(x), schneller als cl_coshsinh_ratseries
				var cl_LF xx = extend(x,TheLfloat(x)->len+1);
				var cl_F y = exp(xx);
				var cl_F z = scale_float(y + recip(y), -1); // (/ (+ y (/ y)) 2)
				return cl_float(z,x);
			} else {
				var cl_LF xx = The(cl_LF)(cl_F_extendsqrt(x));
				var cl_LF y = scale_float(xx,-1);
				// 1 + 2*sinh(y)^2, und wieder runden
				return cl_float(1 + scale_float(sinhx_naive(y),1), x);
			}
		} else {
			var cl_F xx = cl_F_extendsqrt(x);
			var cl_F y = scale_float(xx,-1);
			// 1 + 2*y^2*(sinh(y)/y)^2, und wieder runden
			return cl_float(1 + scale_float(square(y) * sinhxbyx_naive(y),1), x);
		}
	} else {
		// e>=0 -> verwende exp(x)
		var cl_F y = exp(x);
		return scale_float(y + recip(y), -1); // (/ (+ y (/ y)) 2)
	}
}

// Timings of the three algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N      naive  ratseries exp&recip
//   10     0.008   0.037     0.012
//   25     0.032   0.117     0.047
//   50     0.11    0.33      0.017
//  100     0.40    1.06      0.63
//  250     2.65    5.2       3.3
//  500    11.1    18.7      11.5
// 1000    46      61        35
// 2500   238     250       143
// ==> exp&recip fastest for N >= 600.

}  // namespace cln
