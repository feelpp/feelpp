// expx().

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
#include "float/lfloat/misc/cl_LF_exponent.cc"

namespace cln {

// cl_F expx_naive (const cl_F& x)
// cl_LF expx_naive (const cl_LF& x)
//
// Methode:
// e := Exponent aus (decode-float x), d := (float-digits x)
// Bei x=0.0 oder e<-d liefere 1.0
//   (denn bei e<=-d-1 ist abs(exp(x)-1) = abs(x)+O(x^2) < 2^(-d-1),
//    also ist exp(x), auf d Bits gerundet, gleich 1.0).
// Bei e<=-sqrt(d) verwende die Potenzreihe
//   exp(x) = sum(j=0..inf,x^j/j!):
//   b:=1, i:=0, sum:=0,
//   while (/= sum (setq sum (+ sum b))) do b:=b*x/(i+1), i:=i+1.
//   Ergebnis sum.
// Sonst setze y := x/2 = (scale-float x -1),
//   berechne rekursiv z:=exp(y) und liefere z^2.
// Aufwand: asymptotisch d^2.5 .

const cl_LF expx_naive (const cl_LF& x)
{
// Methode:
// wie oben, mit adaptiver Genauigkeit während der Potenzreihen-Summation.
	if (zerop_inline(x))
		return cl_float(1,x);
	var uintC actuallen = TheLfloat(x)->len;
	var uintC d = float_digits(x);
	var sintE e = float_exponent_inline(x);
	if (e < -(sintC)d) // e < -d ?
		return cl_float(1,x); // ja -> 1.0 als Ergebnis
 {	Mutable(cl_LF,x);
	var uintE k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist nicht schlecht,
	// auch im Bereich d = ca. 800.
	var sintL e_limit = -1-isqrtC(d); // -1-floor(sqrt(d))
	if (e > e_limit) {
		// e > -1-floor(sqrt(d)) -> muß |x| verkleinern.
		k = e - e_limit;
		x = scale_float(x,-(sintE)k); // x := x/2^k
		// Neuer Exponent = e-k = e_limit.
	}
	// Potenzreihe anwenden:
	var int i = 0;
	var cl_LF b = cl_float(1,x); // b := (float 1 x)
	var cl_LF eps = scale_float(b,-(sintC)d-10);
	var cl_LF sum = cl_float(0,x); // sum := (float 0 x)
	loop {
		var cl_LF new_sum = sum + LF_to_LF(b,actuallen);
		if (new_sum == sum) // = sum ?
			break; // ja -> Potenzreihe abbrechen
		sum = new_sum;
		i = i+1;
		b = cl_LF_shortenwith(b,eps);
		b = (b*x)/(cl_I)i; // b := b*x/i
	}
	var cl_LF& result = sum; // sum als Ergebnis
	// Wegen Rekursion noch k mal quadrieren:
	for ( ; k > 0; k--)
		result = square(result);
	return result;
}}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

const cl_F expx_naive (const cl_F& x)
{
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		return expx_naive(x);
	}
	if (zerop(x))
		return cl_float(1,x);
	var uintC d = float_digits(x);
	var sintE e = float_exponent(x);
	if (e < -(sintC)d) // e < -d ?
		return cl_float(1,x); // ja -> 1.0 als Ergebnis
 {	Mutable(cl_F,x);
	var uintE k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist nicht schlecht. Für
	// d > 1600 scheint der Bereich 2.0 <= limit_slope <= 2.6 am besten
	// zu sein (mit bis zu 15% Beschleunigung gegenüber limit_slope = 1.0),
	// aber in diesem Bereich rechnen wir gar nicht.
	// Wir wählen limit_slope = 1.5.
	var sintL e_limit = -1-floor(isqrtC(d)*3,2); // -1-floor(sqrt(d))
	if (e > e_limit) {
		// e > -1-floor(sqrt(d)) -> muß |x| verkleinern.
		k = e - e_limit;
		x = scale_float(x,-(sintE)k); // x := x/2^k
		// Neuer Exponent = e-k = e_limit.
	}
	// Potenzreihe anwenden:
	var int i = 0;
	var cl_F b = cl_float(1,x); // b := (float 1 x)
	var cl_F sum = cl_float(0,x); // sum := (float 0 x)
	loop {
		var cl_F new_sum = sum + b;
		if (new_sum == sum) // = sum ?
			break; // ja -> Potenzreihe abbrechen
		sum = new_sum;
		i = i+1;
		b = (b*x)/(cl_I)i; // b := b*x/i
	}
	var cl_F& result = sum; // sum als Ergebnis
	// Wegen Rekursion noch k mal quadrieren:
	for ( ; k > 0; k--)
		result = square(result);
	return result;
}}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

const cl_LF expx_ratseries (const cl_LF& x)
{
	// [Jonathan M. Borwein, Peter B. Borwein: Pi and the AGM.
	//  Wiley 1987. Section 10.2.3]
	var uintC len = TheLfloat(x)->len;
	var cl_idecoded_float x_ = integer_decode_float(x);
	// x = (-1)^sign * 2^exponent * mantissa
	var uintE lq = cl_I_to_UE(- x_.exponent);
	var const cl_I& p = x_.mantissa;
	// Compute exp(p/2^lq) by splitting into pieces.
	// Each piece gives rise to a factor exp(pk/2^lqk).
	// Instead of the standard choice lqk = 2^k, we choose
	// lqk = c^k + O(1), where c > 1 is real.
	// Running time on Linux i486, 33 Mhz, computing exp(sqrt(2)-1):
	//   c  2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5
	//  (a) 400 393 390 377 371 360 363 367 367 358 362 362 363 362 376 372
	//  (b) 311 317 305 312 295 291 286 293 291 284 295 284 293 287 288 305
	// (a): N=300, time in 0.01 sec. (b): N=1000, time in 0.1 sec.
	// Values 2.5 <= c <= 3.2 seem best. Let's choose c = 2.875.
	var bool first_factor = true;
	var cl_LF product;
	var uintE b1;
	var uintE b2;
	for (b1 = 0, b2 = 1; b1 < lq; b1 = b2, b2 = ceiling(b2*23,8)) {
		// Piece containing bits b1+1..b2 after "decimal point"
		// in the binary representation of (p/2^lq).
		var uintE lqk = (lq >= b2 ? b2 : lq);
		var cl_I pk = ldb(p,cl_byte(lqk-b1,lq-lqk));
		// Compute exp(pk/2^lqk).
		if (!zerop(pk)) {
			if (minusp(x_.sign)) { pk = -pk; }
			var cl_LF factor = cl_exp_aux(pk,lqk,len);
			if (first_factor) {
				product = factor;
				first_factor = false;
			} else
				product = product * factor;
		}
	}
	if (first_factor)
		return cl_I_to_LF(1,len);
	else
		return product;
}
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
// ("naive" with adaptive limit_slope, about sqrt(ln(len)).)
//   N      naive  ratseries
//   10     0.010   0.027
//   25     0.039   0.072
//   50     0.15    0.19
//  100     0.60    0.55
//  250     3.9     2.6
//  500    16.3     9.3
// 1000    68      29
// ==> ratseries faster for N >= 84.

}  // namespace cln
