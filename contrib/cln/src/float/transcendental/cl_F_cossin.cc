// cos_sin().

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

const cos_sin_t cos_sin (const cl_F& x)
{
// Methode:
// Genauigkeit erhöhen,
// (q,r) := (round x (float pi/2 x)), so daß |r|<=pi/4.
// y:=(sin(r)/r)^2 errechnen.
// cos(r) berechnen:
//   e := Exponent aus (decode-float r), d := (float-digits r)
//   Bei r=0.0 oder e<=-d/2 liefere 1.0
//     (denn bei e<=-d/2 ist r^2/2 < 2^(-d)/2 = 2^(-d-1), also
//     1 >= cos(r) > 1-r^2/2 > 1-2^(-d-1),
//     also ist cos(r), auf d Bits gerundet, gleich 1.0).
//   Sonst sqrt(1-r^2*y).
// sin(r) berechnen: r*sqrt(y).
// Genauigkeit wieder verringern.
// Falls q = 0 mod 4: (cos(r), sin(r))
// Falls q = 1 mod 4: (-sin(r), cos(r))
// Falls q = 2 mod 4: (-cos(r), -sin(r))
// Falls q = 3 mod 4: (sin(r), -cos(r))

	// Rechengenauigkeit erhöhen und durch pi/2 dividieren:
	var cl_F cos_r;
	var cl_F sin_r;
	var cl_I q;
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		if (TheLfloat(x)->len >= 2710) {
			var cl_F_div_t q_r = cl_round_pi2(extend(x,TheLfloat(x)->len+1));
			q = q_r.quotient;
			var cl_LF r = The(cl_LF)(q_r.remainder);
			var cl_LF_cos_sin_t trig = cl_cossin_ratseries(r);
			cos_r = cl_float(trig.cos,x);
			sin_r = cl_float(trig.sin,x);
		} else {
			var cl_F_div_t q_r = cl_round_pi2(cl_F_extendsqrt(x));
			q = q_r.quotient;
			var cl_LF r = The(cl_LF)(q_r.remainder);
			var cl_LF y = sinx_naive(r); // y := sin(r)^2
			// erste Komponente cos(r) berechnen:
			if (zerop(r) || (float_exponent(r) <= (-(sintC)float_digits(r))>>1))
				cos_r = cl_float(1,x); // cos(r) = 1.0
			else
				cos_r = cl_float(sqrt(1-y),x); // cos(r) = sqrt(1-y)
			// zweite Komponente sin(r) berechnen:
			sin_r = cl_float(sqrt(y),x);
			if (minusp(r))
				sin_r = - sin_r;
		}
	} else {
		var cl_F_div_t q_r = cl_round_pi2(cl_F_extendsqrt(x));
		q = q_r.quotient;
		var cl_F& r = q_r.remainder;
		var cl_F y = sinxbyx_naive(r); // y := (sin(r)/r)^2
		// erste Komponente cos(r) berechnen:
		if (zerop(r) || (float_exponent(r) <= (-(sintC)float_digits(r))>>1))
			cos_r = cl_float(1,x); // cos(r) = 1.0
		else
			cos_r = cl_float(sqrt(1 - square(r)*y),x); // sqrt(1-r^2*y)
		// zweite Komponente sin(r) berechnen:
		sin_r = cl_float(r*sqrt(y),x);
	}
	// evtl. Vorzeichenwechsel oder Vertauschen:
	switch (cl_I_to_UL(logand(q,3))) { // q mod 4
		case 0: return cos_sin_t(cos_r,sin_r);
		case 1: return cos_sin_t(-sin_r,cos_r);
		case 2: return cos_sin_t(-cos_r,-sin_r);
		case 3: return cos_sin_t(sin_r,-cos_r);
		default: NOTREACHED
	}
}

}  // namespace cln
