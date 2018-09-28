// sin().

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

const cl_F sin (const cl_F& x)
{
// Methode:
// Genauigkeit erhöhen,
// (q,r) := (round x (float pi/2 x)), so daß |r|<=pi/4.
// y:=(sin(r)/r)^2 errechnen.
// Falls q gerade:
//   sin(r) berechnen: r*sqrt(y).
// Falls q ungerade:
//   cos(r) berechnen:
//     e := Exponent aus (decode-float r), d := (float-digits r)
//     Bei r=0.0 oder e<=-d/2 liefere 1.0
//       (denn bei e<=-d/2 ist r^2/2 < 2^(-d)/2 = 2^(-d-1), also
//       1 >= cos(r) > 1-r^2/2 > 1-2^(-d-1),
//       also ist cos(r), auf d Bits gerundet, gleich 1.0).
//     Sonst sqrt(1-r^2*y).
// Falls q == 2,3 mod 4, Vorzeichenwechsel.

	// Rechengenauigkeit erhöhen und durch pi/2 dividieren:
	var cl_F z;
	var cl_I q;
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		if (TheLfloat(x)->len >= 2750) {
			var cl_F_div_t q_r = cl_round_pi2(extend(x,TheLfloat(x)->len+1));
			q = q_r.quotient;
			var cl_LF r = The(cl_LF)(q_r.remainder);
			var cl_LF_cos_sin_t trig = cl_cossin_ratseries(r);
			if (evenp(q))
				z = cl_float(trig.sin,x);
			else
				z = cl_float(trig.cos,x);
		} else {
			var cl_F_div_t q_r = cl_round_pi2(cl_F_extendsqrt(x));
			q = q_r.quotient;
			var cl_LF r = The(cl_LF)(q_r.remainder);
			var cl_LF y = sinx_naive(r); // y := sin(r)^2
			if (evenp(q)) {
				// sin(r) berechnen:
				z = cl_float(sqrt(y),x);
				if (minusp(r))
					z = -z;
			} else {
				// cos(r) berechnen:
				if (zerop(r) || (float_exponent(r) <= (-(sintC)float_digits(r))>>1))
					z = cl_float(1,x); // cos(r) = 1.0
				else
					z = cl_float(sqrt(1 - y),x); // sqrt(1-y)
			}
		}
	} else {
		var cl_F_div_t q_r = cl_round_pi2(cl_F_extendsqrt(x));
		q = q_r.quotient;
		var cl_F& r = q_r.remainder;
		var cl_F y = sinxbyx_naive(r); // y := (sin(r)/r)^2
		if (evenp(q)) {
			// sin(r) berechnen:
			z = cl_float(r*sqrt(y),x);
		} else {
			// cos(r) berechnen:
			if (zerop(r) || (float_exponent(r) <= (-(sintC)float_digits(r))>>1))
				z = cl_float(1,x); // cos(r) = 1.0
			else
				z = cl_float(sqrt(1 - square(r)*y),x); // sqrt(1-r^2*y)
		}
	}
	// evtl. Vorzeichenwechsel:
	if (cl_I_to_UL(logand(q,2))==0)
		return z;
	else
		return -z;
}

// Timings of the two algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N      naive  ratseries
//   10     0.010   0.048
//   25     0.035   0.119
//   50     0.12    0.37
//  100     0.44    1.09
//  250     2.8     5.5
//  500    11.6    19.4
// 1000    48      64
// 2500   243     261
// ==> ratseries faster for N >= 2750.

}  // namespace cln
