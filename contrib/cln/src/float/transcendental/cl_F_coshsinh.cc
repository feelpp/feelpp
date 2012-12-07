// cosh_sinh().

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

const cosh_sinh_t cosh_sinh (const cl_F& x)
{
// Methode:
// Genauigkeit erhöhen,
// e := Exponent aus (decode-float x), d := (float-digits x)
// falls x=0.0 oder e<=(1-d)/2 liefere (1.0,x)
//   (denn bei e<=(1-d)/2 ist
//    1 <= sinh(x)/x < cosh(x) = 1+x^2/2+... < 1+2^(-d),
//    also ist cosh(x), auf d Bits gerundet, gleich 1.0
//    und sinh(x), auf d Bits gerundet, gleich x).
// falls e<0:
//   y:=(sinh(x)/x)^2 errechnen,
//   cosh(x) = sqrt(1+x^2*y) und sinh(x) = x*sqrt(y) errechnen.
// falls e>=0: y:=exp(x) errechnen,
//   (scale-float (+ y (/ y)) -1) und (scale-float (- y (/ y)) -1) bilden.
// Genauigkeit wieder verringern.

	var sintE e = float_exponent(x);
	if (e < 0) { // Exponent e abtesten
		// e<0
		if (zerop(x) || (e <= (1-(sintC)float_digits(x))>>1))
			// e <= (1-d)/2 <==> e <= -ceiling((d-1)/2)
			return cosh_sinh_t(cl_float(1,x),x);
		// Rechengenauigkeit erhöhen
		if (longfloatp(x)) {
			DeclareType(cl_LF,x);
			#if 0
			if (TheLfloat(x)->len >= infty) {
				var cl_LF xx = extend(x,TheLfloat(x)->len+1);
				var cl_LF_cosh_sinh_t hyp = cl_coshsinh_ratseries(xx);
				return cosh_sinh_t(
					cl_float(hyp.cosh,x),
					cl_float(hyp.sinh,x)
				       );
			} else
			#endif
			if (TheLfloat(x)->len >= 585) {
				// verwende exp(x), schneller als cl_coshsinh_ratseries
				var cl_LF xx = extend(x,TheLfloat(x)->len+ceiling((uintE)(-e),intDsize));
				var cl_F y = exp(xx);
				var cl_F y_inv = recip(y);
				return cosh_sinh_t(
					cl_float(scale_float(y + y_inv, -1), x),
					cl_float(scale_float(y - y_inv, -1), x)
				       );
			} else {
				var cl_LF xx = The(cl_LF)(cl_F_extendsqrt(x));
				var cl_LF y = sinhx_naive(xx);
				var cl_LF z = sqrt(y);
				if (minusp(xx))
					z = -z;
				return cosh_sinh_t(
					cl_float(sqrt(1+y),x), // sqrt(1+y)
					cl_float(z,x)
				       );
			}
		} else {
			var cl_F xx = cl_F_extendsqrt(x);
			var cl_F y = sinhxbyx_naive(xx);
			return cosh_sinh_t(
				cl_float(sqrt(1+square(xx)*y),x), // sqrt(1+x^2*y)
				cl_float(xx*sqrt(y),x)
			       );
		}
	} else {
		// e>=0 -> verwende exp(x)
		var cl_F y = exp(x);
		var cl_F y_inv = recip(y);
		return cosh_sinh_t(
			scale_float(y+y_inv,-1),
			scale_float(y-y_inv,-1)
		       );
	}
}

}  // namespace cln
