// atanx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/float.h"
#include "base/cl_low.h"
#include "float/cl_F.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"
#include "cln/integer.h"

#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_zerop.cc"
#include "float/lfloat/elem/cl_LF_minusp.cc"
#include "float/lfloat/misc/cl_LF_exponent.cc"

namespace cln {

// cl_F atanx_naive (const cl_F& x)
// cl_LF atanx_naive (const cl_LF& x)
//
// Methode:
// e := Exponent aus (decode-float x), d := (float-digits x)
// Bei x=0.0 oder e<=-d/2 liefere x
//   (denn bei e<=-d/2 ist x^2/3 < x^2/2 < 2^(-d)/2 = 2^(-d-1), also
//   1 >= atan(x)/x > 1-x^2/3 > 1-2^(-d-1),
//   also ist atan(x)/x, auf d Bits gerundet, gleich 1.0).
// Bei e<=-sqrt(d) verwende die Potenzreihe
//   atan(x)/x = sum(j=0..inf,(-x^2)^j/(2j+1)):
//   a:=-x^2, b:=1, i:=1, sum:=0,
//   while (/= sum (setq sum (+ sum (/ b i)))) do i:=i+2, b:=b*a.
//   Ergebnis x*sum.
// Sonst setze y := x/(1+sqrt(1+x^2)), berechne rekursiv z:=atan(y)
//   und liefere 2*z = (scale-float z 1).
// Diese Rekursion wird entrekursiviert. Statt k mal hintereinander
//   x := x/(1+sqrt(1+x^2)) zu bilden, arbeitet man lieber mit den Kehrwerten,
//   setzt also x := 1/|x|, dann k mal x := x+sqrt(x^2+1), dann x := +- 1/x.
// Aufwand: asymptotisch d^2.5 .

static const cl_LF atanx_naive (const cl_LF& x)
{
	if (zerop_inline(x))
		return x;
	var uintC actuallen = TheLfloat(x)->len;
	var uintC d = float_digits(x);
	var sintE e = float_exponent_inline(x);
	if (e <= (sintC)(-d)>>1) // e <= -d/2 <==> e <= -ceiling(d/2)
		return x; // ja -> x als Ergebnis
	var uintL k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist schlecht (ca. 20% zu
	// schlecht). Ein guter Wert ist:
	// Für naive1: limit_scope = 0.5.
	// Für naive2: limit_scope = 0.375 (ca. 0.5 für kleine len, 0.35 für
	// große len).
	var uintL sqrt_d = floor(isqrtC(d)*3,8); // limit_slope*floor(sqrt(d))
	var cl_LF xx = x;
	if (e >= (sintL)(-sqrt_d)) {
		// e > -1-limit_slope*floor(sqrt(d)) -> muß |x| verkleinern.
		var sintL e_limit = 1+sqrt_d; // 1+limit_slope*floor(sqrt(d))
		xx = recip(abs(xx)); // 1/|x|
		do {
		  // nächstes x nach der Formel x := x+sqrt(x^2 + 1) berechnen:
		  xx = sqrt(square(xx) + cl_float(1,xx)) + xx;
		  k = k+1;
		} until (float_exponent_inline(xx) > e_limit);
		// Schleifenende mit Exponent(x) > 1+limit_slope*floor(sqrt(d)),
		// also x >= 2^(1+limit_slope*floor(sqrt(d))),
		// also 1/x <= 2^(-1-limit_slope*floor(sqrt(d))).
		// Nun kann die Potenzreihe auf 1/x angewandt werden.
		xx = recip(xx);
		if (minusp_inline(x))
			xx = - xx; // Vorzeichen wieder rein
	}
	// Potenzreihe anwenden:
	var int i = 1;
	var cl_LF a = - square(xx); // a = - x^2
	var cl_LF b = cl_float(1,xx); // b := (float 1 x)
	var cl_LF sum = cl_float(0,xx); // sum := (float 0 x)
	if (0) {
		// naive1:
		// floating-point representation
		loop {
			var cl_LF new_sum = sum + b / (cl_I)i; // (+ sum (/ b i))
			if (new_sum == sum) // = sum ?
				break; // ja -> Potenzreihe abbrechen
			sum = new_sum;
			b = b*a;
			i = i+2;
		}
	} else {
		// naive2:
		// floating-point representation with smooth precision reduction
		var cl_LF eps = scale_float(b,-(sintC)d-10);
		loop {
			var cl_LF new_sum = sum + LF_to_LF(b/(cl_I)i,actuallen); // (+ sum (/ b i))
			if (new_sum == sum) // = sum ?
				break; // ja -> Potenzreihe abbrechen
			sum = new_sum;
			b = cl_LF_shortenwith(b,eps);
			b = b*a;
			i = i+2;
		}
	}
	var cl_LF erg = sum*xx; // sum*x als Ergebnis
	return scale_float(erg,k); // wegen Rekursion noch mal 2^k
}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

static const cl_F atanx_naive (const cl_F& x)
{
	if (zerop(x))
		return x;
	var uintC d = float_digits(x);
	var sintE e = float_exponent(x);
	if (e <= (sintC)(-d)>>1) // e <= -d/2 <==> e <= -ceiling(d/2)
		return x; // ja -> x als Ergebnis
	var uintL k = 0; // Rekursionszähler k:=0
	var uintL sqrt_d = floor(isqrtC(d),2); // limit_slope*floor(sqrt(d))
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist schlecht (ca. 20% zu
	// schlecht). Ein guter Wert ist limit_scope = 0.5.
	var cl_F xx = x;
	if (e >= (sintL)(-sqrt_d)) {
		// e > -1-limit_slope*floor(sqrt(d)) -> muß |x| verkleinern.
		var sintL e_limit = 1+sqrt_d; // 1+limit_slope*floor(sqrt(d))
		xx = recip(abs(xx)); // 1/|x|
		do {
		  // nächstes x nach der Formel x := x+sqrt(x^2 + 1) berechnen:
		  xx = sqrt(square(xx) + cl_float(1,xx)) + xx;
		  k = k+1;
		} until (float_exponent(xx) > e_limit);
		// Schleifenende mit Exponent(x) > 1+limit_slope*floor(sqrt(d)),
		// also x >= 2^(1+limit_slope*floor(sqrt(d))),
		// also 1/x <= 2^(-1-limit_slope*floor(sqrt(d))).
		// Nun kann die Potenzreihe auf 1/x angewandt werden.
		xx = recip(xx);
		if (minusp(x))
			xx = - xx; // Vorzeichen wieder rein
	}
	// Potenzreihe anwenden:
	var int i = 1;
	var cl_F a = - square(xx); // a = - x^2
	var cl_F b = cl_float(1,xx); // b := (float 1 x)
	var cl_F sum = cl_float(0,xx); // sum := (float 0 x)
	loop {
		var cl_F new_sum = sum + b / (cl_I)i; // (+ sum (/ b i))
		if (new_sum == sum) // = sum ?
			break; // ja -> Potenzreihe abbrechen
		sum = new_sum;
		b = b*a;
		i = i+2;
	}
	var cl_F erg = sum*xx; // sum*x als Ergebnis
	return scale_float(erg,k); // wegen Rekursion noch mal 2^k
}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

static const cl_LF atanx_ratseries (const cl_LF& t)
{
	// Method:
	// Based on the same ideas as lnx_ratseries.
	//   e := exponent of (decode-float t), d := (float-digits t).
	//   If t=0.0 or e<=-d/2, return t.
	// (x,y) := (1/sqrt(1+t^2),t/sqrt(1+t^2)), z := 0.
	// Loop
	//   [(x+i*y)*exp(i*z) is invariant, x>0, sqrt(x^2+y^2)=1]
	//   e := exponent of (decode-float y), d := (float-digits y).
	//   If y=0.0 or e<=-d/2, return z+y
	//   (because if e<=-d/2 then |y|^3/6 < 2^(-d)/2*|y|, and since
	//   asin(y) = y+y^3/6+..., asin(y) rounded to d bits is = y).
	//   Choose approximation z' of angle(x+i*y):
	//     If |y| >= 1/2, set z' = 1/2 * sign(y).
	//     If |y| < 2^-n with n maximal, set
	//       z' = truncate(y*2^(2n))/2^(2n).
	//   Set z := z + z' and x+i*y := (x+i*y)*exp(-i*z').
	var uintC len = TheLfloat(t)->len;
	var uintC d = intDsize*len;
	if (zerop_inline(t) || (float_exponent_inline(t) <= (sintC)(-d)>>1))
		return t;
	var cl_LF x = recip(sqrt(cl_I_to_LF(1,len) + square(t)));
	var cl_LF y = t*x;
	var cl_LF z = cl_I_to_LF(0,len);
	loop {
		if (zerop_inline(y) || (float_exponent_inline(y) <= (sintC)(-d)>>1))
			break;
		var cl_idecoded_float y_ = integer_decode_float(y);
		// y = (-1)^sign * 2^exponent * mantissa
		var uintC lm = integer_length(y_.mantissa);
		var uintE me = cl_I_to_UE(- y_.exponent);
		var cl_I p;
		var uintE lq;
		var bool last_step = false;
		if (lm >= me) { // |y| >= 1/2 ?
			p = y_.sign; // 1 or -1
			lq = 1;
		} else {
			var uintE n = me - lm; // |y| < 2^-n with n maximal
			// Set p to the first n bits of |y|:
			if (lm > n) {
				p = y_.mantissa >> (lm - n);
				lq = 2*n;
			} else {
				p = y_.mantissa;
				lq = lm + n;
			}
			if (minusp(y_.sign)) { p = -p; }
			// If 2*n >= lm = intDsize*len, then within our
			// precision exp(-i*z')=1-i*z' (because |z'^2| < 2^-lm),
			// and we know a priori that the iteration will stop
			// after the next big multiplication. This saves one
			// big multiplication at the end.
			if (2*n >= lm)
				last_step = true;
		}
		z = z + scale_float(cl_I_to_LF(p,len),-(sintE)lq);
		if (last_step)
			break;
		var cl_LF_cos_sin_t cis_z = cl_cossin_aux(-p,lq,len);
		var cl_LF new_x = x*cis_z.cos - y*cis_z.sin;
		var cl_LF new_y = x*cis_z.sin + y*cis_z.cos;
		x = new_x; y = new_y;
	}
	return z+y;
}
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N      naive1  naive2  ratseries
//   10     0.013   0.013   0.043
//   25     0.062   0.048   0.122
//   50     0.25    0.17    0.34
//  100     1.06    0.70    1.07
//  250     7.5     5.0     5.6
//  500    34.7    23.2    20.0
// 1000   167     112      65
// ==> ratseries faster for N >= 325.

const cl_F CL_FLATTEN atanx (const cl_F& x)
{
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		if (TheLfloat(x)->len >= 325)
			return cl_float(atanx_ratseries(extend(x,TheLfloat(x)->len+1)),x);
		else
			return atanx_naive(x);
	} else
		return atanx_naive(x);
}
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

}  // namespace cln
