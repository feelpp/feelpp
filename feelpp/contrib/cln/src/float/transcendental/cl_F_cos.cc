// cos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/transcendental/cl_F_tran.h"
#include "float/cl_F.h"
#include "cln/integer.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F cos (const cl_F& x)
{
// Methode:
// Genauigkeit erhöhen,
// (q,r) := (round x (float pi x)), so daß |r|<=pi/2.
// e := Exponent aus (decode-float r), d := (float-digits r)
// Bei r=0.0 oder e<=-d/2 liefere 1.0
//   (denn bei e<=-d/2 ist r^2/2 < 2^(-d)/2 = 2^(-d-1), also
//   1 >= cos(r) > 1-r^2/2 > 1-2^(-d-1),
//   also ist cos(r), auf d Bits gerundet, gleich 1.0).
// Sonst s := r/2 = (scale-float r -1),
//   (sin(s)/s)^2 errechnen, cos(r) = 1-r*s*(sin(s)/s)^2 errechnen.
// Falls q ungerade: Vorzeichenwechsel.

	// Rechengenauigkeit erhöhen und durch pi dividieren:
	var cl_F cos_r;
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		if (TheLfloat(x)->len >= 2850) {
			var cl_F_div_t q_r = cl_round_pi2(extend(x,TheLfloat(x)->len+1));
			var cl_I& q = q_r.quotient;
			var cl_LF r = The(cl_LF)(q_r.remainder);
			var cl_LF_cos_sin_t trig = cl_cossin_ratseries(r);
			switch (cl_I_to_UL(logand(q,3))) { // q mod 4
				case 0: return cl_float(trig.cos,x);
				case 1: return -cl_float(trig.sin,x);
				case 2: return -cl_float(trig.cos,x);
				case 3: return cl_float(trig.sin,x);
				default: NOTREACHED
			}
		} else {
			var cl_F_div_t q_r = cl_round_pi(cl_F_extendsqrt(x));
			var cl_I& q = q_r.quotient;
			var cl_LF r = The(cl_LF)(q_r.remainder);
			if (zerop(r) || (float_exponent(r) <= (-(sintC)float_digits(r))>>1))
				cos_r = cl_float(1,x); // (cos r) = 1.0
			else {
				var cl_LF s = scale_float(r,-1); // s := r/2
				cos_r = cl_float(1-scale_float(sinx_naive(s),1),x); // cos(2s) = 1-2*sin(s)^2
			}
			if (oddp(q))
				return -cos_r; // q ungerade -> mal -1
			else
				return cos_r;
		}
	} else {
		var cl_F_div_t q_r = cl_round_pi(cl_F_extendsqrt(x));
		var cl_I& q = q_r.quotient;
		var cl_F& r = q_r.remainder;
		if (zerop(r) || (float_exponent(r) <= (-(sintC)float_digits(r))>>1))
			cos_r = cl_float(1,x); // (cos r) = 1.0
		else {
			var cl_F s = scale_float(r,-1); // s := r/2
			cos_r = cl_float(1 - r * s * sinxbyx_naive(s),x);
		}
		if (oddp(q))
			return -cos_r; // q ungerade -> mal -1
		else
			return cos_r;
	}
}

// Timings of the two algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N      naive  ratseries
//   10     0.009   0.049
//   25     0.033   0.137
//   50     0.11    0.37
//  100     0.41    1.15
//  250     2.7     5.5
//  500    11.1    19.4
// 1000    46      64
// 2500   239     260
// ==> ratseries faster for N >= 2850.

}  // namespace cln
