// atanhx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"
#include "float/cl_F.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/float.h"
#include "base/cl_low.h"

#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_zerop.cc"
#include "float/lfloat/elem/cl_LF_minusp.cc"
#include "float/lfloat/misc/cl_LF_exponent.cc"

namespace cln {

// cl_F atanhx (const cl_F& x)
// cl_LF atanhx (const cl_LF& x)
//
// Methode:
// e := Exponent aus (decode-float x), d := (float-digits x)
// Bei x=0.0 oder e<=-d/2 liefere x
//   (denn bei e<=-d/2 ist x^2 < 2^(-d), also
//   1 <= atanh(x)/x = 1+x^2/3+x^4/5+... < 1+x^2/2 < 1+2^(-d-1) < 1+2^(-d),
//   also ist atanh(x)/x, auf d Bits gerundet, gleich 1.0).
// Bei großem d verwende die Formel ln((1+x)/(1-x))/2 (asymptotisch schneller),
//   aber erhöhe die Genauigkeit, so daß beim Bilden von 1+x keine Bits verloren
//   gehen.
// Bei e<=-sqrt(d) verwende die Potenzreihe
//   atanh(x)/x = sum(j=0..inf,(x^2)^j/(2j+1)):
//   a:=x^2, b:=1, i:=1, sum:=0,
//   while (/= sum (setq sum (+ sum (/ b i)))) do i:=i+2, b:=b*a.
//   Ergebnis x*sum.
// Sonst setze y := x/(1+sqrt(1-x^2)), berechne rekursiv z:=atanh(y)
//   und liefere 2*z = (scale-float z 1).
// Diese Rekursion wird entrekursiviert. Statt k mal hintereinander
//   x := x/(1+sqrt(1-x^2)) zu bilden, arbeitet man lieber mit den Kehrwerten,
//   setzt also x := 1/|x|, dann k mal x := x+sqrt(x^2-1), dann x := +- 1/x.
// Aufwand: asymptotisch d^2.5 .

const cl_LF atanhx (const cl_LF& x)
{
	if (zerop_inline(x))
		return x;
	var uintC actuallen = TheLfloat(x)->len;
	var uintC d = float_digits(x);
	var sintE e = float_exponent_inline(x);
	if (e <= (sintC)(-d)>>1) // e <= -d/2 <==> e <= -ceiling(d/2)
		return x; // ja -> x als Ergebnis
	if (actuallen >= 34) {
		DeclareType(cl_LF,x);
		var cl_LF xx = extend(x,TheLfloat(x)->len+ceiling((uintE)(-e),intDsize));
		return cl_float(scale_float(ln((1+xx)/(1-xx)),-1),x);
	}
	var uintL k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist schlecht (ca. 15% zu
	// schlecht). Ein guter Wert ist:
	// für naive1: limit_scope = 0.625 = 5/8,
	// für naive2: limit_scope = 0.4 = 13/32.
	var uintL sqrt_d = floor(isqrtC(d)*13,32); // limit_slope*floor(sqrt(d))
	var cl_LF xx = x;
	if (e >= (sintL)(-sqrt_d)) {
		// e > -1-limit_slope*floor(sqrt(d)) -> muß |x| verkleinern.
		var sintL e_limit = 1+sqrt_d; // 1+limit_slope*floor(sqrt(d))
		xx = recip(abs(xx)); // 1/|x|
		do {
		  // nächstes x nach der Formel x := x+sqrt(x^2 - 1) berechnen:
		  xx = sqrt(square(xx) + cl_float(-1,xx)) + xx;
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
	var cl_LF a = square(xx); // a = x^2
	var cl_LF b = cl_float(1,xx); // b := (float 1 x)
	var cl_LF sum = cl_float(0,xx); // sum := (float 0 x)
	if (0) {
		// naive1:
		// floating-point representation
		loop {
			var cl_LF new_sum = sum + b/(cl_I)i; // (+ sum (/ b i))
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
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

const cl_F atanhx (const cl_F& x)
{
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		return atanhx(x);
	}
	if (zerop(x))
		return x;
	var uintC d = float_digits(x);
	var sintE e = float_exponent(x);
	if (e <= (sintC)(-d)>>1) // e <= -d/2 <==> e <= -ceiling(d/2)
		return x; // ja -> x als Ergebnis
	var uintL k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist schlecht (ca. 15% zu
	// schlecht). Ein guter Wert ist limit_scope = 0.625 = 5/8.
	var uintL sqrt_d = floor(isqrtC(d)*5,8); // limit_slope*floor(sqrt(d))
	var cl_F xx = x;
	if (e >= (sintL)(-sqrt_d)) {
		// e > -1-limit_slope*floor(sqrt(d)) -> muß |x| verkleinern.
		var sintL e_limit = 1+sqrt_d; // 1+limit_slope*floor(sqrt(d))
		xx = recip(abs(xx)); // 1/|x|
		do {
		  // nächstes x nach der Formel x := x+sqrt(x^2 - 1) berechnen:
		  xx = sqrt(square(xx) + cl_float(-1,xx)) + xx;
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
	var cl_F a = square(xx); // a = x^2
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
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N      naive1  naive2  use ln
//   10     0.013   0.013   0.015
//   25     0.064   0.050   0.049
//   50     0.25    0.018   0.17
//  100     1.07    0.75    0.64
//  250     7.6     5.2     2.7
//  500    35.5    24.2     9.7
// 1000   168     116      29.6
// ==> using ln faster for N >= 34.

}  // namespace cln
