// sinh().

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

const cl_F sinh (const cl_F& x)
{
// Methode:
// Genauigkeit erhöhen,
// e := Exponent aus (decode-float x)
// falls e<0: (sinh(x)/x)^2 errechnen, Wurzel ziehen, mit x multiplizieren.
// falls e>=0: y:=exp(x) errechnen, (scale-float (- y (/ y)) -1) bilden.

	if (float_exponent(x) < 0) { // Exponent e abtesten
		// e<0
		// Rechengenauigkeit erhöhen
		if (longfloatp(x)) {
			DeclareType(cl_LF,x);
			#if 0
			if (TheLfloat(x)->len >= infty) {
				var cl_LF xx = extend(x,TheLfloat(x)->len+1);
				var cl_LF_cosh_sinh_t hyp = cl_coshsinh_ratseries(xx);
				return cl_float(hyp.sinh,x);
			} else
			#endif
			if ((TheLfloat(x)->len >= 500)
			    && (float_exponent(x) > (-(sintC)float_digits(x))>>1)) {
				// verwende exp(x), schneller als cl_coshsinh_ratseries
				// (aber nur bei 0 > e > -d/2, denn wir müssen, um
				// Auslöschung zu verhindern, |e| Bits dazunehmen)
				var cl_LF xx = extend(x,TheLfloat(x)->len+ceiling((uintE)(-float_exponent(x)),intDsize));
				var cl_F y = exp(xx);
				var cl_F z = scale_float(y - recip(y), -1); // (/ (- y (/ y)) 2)
				return cl_float(z,x);
			} else {
				var cl_LF xx = The(cl_LF)(cl_F_extendsqrt(x));
				// Wurzel aus sinh(x)^2 bilden
				var cl_LF z = sqrt(sinhx_naive(xx));
				if (minusp(xx))
					z = -z;
				return cl_float(z,x);
			}
		} else {
			var cl_F xx = cl_F_extendsqrt(x);
			// Wurzel aus (sinh(x)/x)^2 mit x multiplizieren und wieder runden
			return cl_float(sqrt(sinhxbyx_naive(xx))*xx,x);
		}
	} else {
		// e>=0 -> verwende exp(x)
		var cl_F y = exp(x);
		return scale_float(y - recip(y), -1); // (/ (- y (/ y)) 2)
	}
}

// Timings of the two algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N      naive  ratseries
//   10     0.008   0.037
//   25     0.034   0.115
//   50     0.13    0.33
//  100     0.50    1.07
//  250     3.3     5.2
//  500    14.2    18.8
// 1000    59      61
// 2500   297     247
// ==> ratseries faster for N >= 1300.

}  // namespace cln
