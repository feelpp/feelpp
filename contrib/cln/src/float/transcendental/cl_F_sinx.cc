// sinxbyx(), sinx().

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

// sinxbyx is mainly for cl_SF, cl_FF, cl_DF, where we want to avoid underflow.

const cl_F sinxbyx_naive (const cl_F& x)
{
// Methode:
// e := Exponent aus (decode-float x), d := (float-digits x)
// Bei x=0.0 oder e<=-d/2 liefere 1.0
//   (denn bei e<=-d/2 ist x^2/6 < x^2/4 < 2^(-d)/4 = 2^(-d-2), also
//   1 >= sin(x)/x > 1-x^2/6 > 1-2^(-d-2), also 1 >= (sin(x)/x)^2 > 1-2^(-d-1),
//   also ist (sin(x)/x)^2, auf d Bits gerundet, gleich 1.0).
// Bei e<=-sqrt(d) verwende die Potenzreihe
//   sin(x)/x = sum(j=0..inf,(-x^2)^j/(2j+1)!):
//   a:=-x^2, b:=1, i:=1, sum:=0,
//   while (/= sum (setq sum (+ sum b))) do b:=b*a/((i+1)*(i+2)), i:=i+2.
//   Ergebnis sum^2.
// Sonst setze y := x/2 = (scale-float x -1),
//   berechne rekursiv z:=(sin(y)/y)^2 und liefere z*(1-y^2*z).
// [Die Grenze sqrt(d) ergibt sich so:
//  Man braucht bei der Potenzreihe mit x=2^-k etwa j Glieder, mit
//  k*j*ln 2 + j*(ln j - 1) = d, und der Aufwand beträgt etwa 2.8*(j/2)
//  Multiplikationen von d-Bit-Zahlen. Bei Halbierungen bis x=2^-k ist der
//  Gesamtaufwand etwa 2*(k+e)+1.4*j(k). Dieses minimieren nach k: Soll sein
//  -1.4 = d/dk j(k) = (d/dj k(j))^-1 = - j^2/(d+j)*ln 2, also j^2=2(d+j),
//  grob j=sqrt(2d) und damit k=sqrt(d).]
// Aufwand: asymptotisch d^2.5 .

	if (zerop(x))
		return cl_float(1,x);
	var uintC d = float_digits(x);
	var sintE e = float_exponent(x);
	if (e <= (-(sintC)d)>>1) // e <= (-d)/2 <==> e <= -ceiling(d/2) ?
		return cl_float(1,x); // ja -> 1.0 als Ergebnis
 {	Mutable(cl_F,x);
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist sehr schlecht (ca. 30%
	// zu schlecht). Gute Werte bei N limbs:
	//    N     limit_slope
	//    5       0.15-0.30
	//   10          0.25
	//   25       0.25-0.35
	//   50       0.35-0.40
	//  100       0.40-0.45
	//  200       0.40-0.45
	// Wähle limit_slope = 13/32 = 0.4.
	var sintL e_limit = -1-floor(isqrtC(d)*13,32); // -1-floor(sqrt(d))
	if (e > e_limit) {
		// e > -1-limit_slope*floor(sqrt(d)) -> muß |x| verkleinern.
		x = scale_float(x,e_limit-e);
		// Neuer Exponent = e_limit.
	}
	var cl_F x2 = square(x);	// x^2
	// Potenzreihe anwenden:
	var cl_F a = - x2; // a := -x^2
	var int i = 1;
	var cl_F b = cl_float(1,x); // b := (float 1 x)
	var cl_F sum = cl_float(0,x); // sum := (float 0 x)
	loop {
		var cl_F new_sum = sum + b;
		if (new_sum == sum) // = sum ?
			break; // ja -> Potenzreihe abbrechen
		sum = new_sum;
		b = (b*a)/(cl_I)((i+1)*(i+2));
		i = i+2;
	}
	var cl_F z = square(sum); // sum^2 als Ergebnis
	while (e > e_limit) {
		z = z - x2 * square(z);
		x2 = scale_float(x2,2); // x^2 := x^2*4
		e--;
	}
	return z;
}}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

const cl_LF sinx_naive (const cl_LF& x)
{
// Methode:
// e := Exponent aus (decode-float x), d := (float-digits x)
// Bei x=0.0 oder e<=-d/2 liefere x
//   (denn bei e<=-d/2 ist x^2/6 < x^2/4 < 2^(-d)/4 = 2^(-d-2), also
//   1 >= sin(x)/x > 1-x^2/6 > 1-2^(-d-2), also ist sin(x)^2, auf d Bits
//   gerundet, gleich x).
// Bei e<=-sqrt(d) verwende die Potenzreihe
//   sin(x) = sum(j=0..inf,x*(-x^2)^j/(2j+1)!):
//   a:=-x^2, b:=x, i:=1, sum:=0,
//   while (/= sum (setq sum (+ sum b))) do b:=b*a/((i+1)*(i+2)), i:=i+2.
//   Ergebnis sum^2.
// Sonst setze y := x/2 = (scale-float x -1),
//   berechne rekursiv z:=sin(y)^2 und liefere 4*z*(1-z) = 1-(1-2*z)^2.
// [Die Grenze sqrt(d) ergibt sich so:
//  Man braucht bei der Potenzreihe mit x=2^-k etwa j Glieder, mit
//  k*j*ln 2 + j*(ln j - 1) = d, und der Aufwand beträgt etwa 2.8*(j/2)
//  Multiplikationen von d-Bit-Zahlen. Bei Halbierungen bis x=2^-k ist der
//  Gesamtaufwand etwa 2*(k+e)+1.4*j(k). Dieses minimieren nach k: Soll sein
//  -1.4 = d/dk j(k) = (d/dj k(j))^-1 = - j^2/(d+j)*ln 2, also j^2=2(d+j),
//  grob j=sqrt(2d) und damit k=sqrt(d).]
// Aufwand: asymptotisch d^2.5 .

	if (zerop_inline(x))
		return x;
	var uintC actuallen = TheLfloat(x)->len;
	var uintC d = float_digits(x);
	var sintE e = float_exponent_inline(x);
	if (e <= (-(sintC)d)>>1) // e <= (-d)/2 <==> e <= -ceiling(d/2) ?
		return square(x); // ja -> x^2 als Ergebnis
 {	Mutable(cl_LF,x);
	var sintE ee = e;
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden. limit_slope = 1.0 ist schlecht (ca. 10% zu
	// schlecht). Ein guter Wert für naive1 ist limit_slope = 0.6,
	// für naive3 aber limit_slope = 0.5.
	var sintL e_limit = -1-floor(isqrtC(d),2); // -1-floor(sqrt(d))
	if (e > e_limit) {
		// e > -1-limit_slope*floor(sqrt(d)) -> muß |x| verkleinern.
		x = scale_float(x,e_limit-e);
		ee = e_limit; // Neuer Exponent = e_limit.
	}
	var cl_LF x2 = square(x); // x^2
	// Potenzreihe anwenden:
	var cl_LF powser_value;
	var cl_LF a = - x2; // a := -x^2
	var int i = 1;
	if (0) {
		// naive1:
		// fixed-point representation
		d = d-ee; // fixed-point representation with d mantissa bits
		var cl_I b = round1(scale_float(x,d)); // b := x
		var cl_I sum = 0; // sum := (float 0 x)
		loop {
			if (b == 0) break;
			sum = sum + b;
			b = round1(round1(The(cl_LF)(b*a)),(cl_I)((i+1)*(i+2)));
			i = i+2;
		}
		powser_value = scale_float(cl_float(sum,x),-(sintC)d);
	} else if (actuallen <= 7) { // Break-even-Point before extendsqrt: N<=6
		// naive2:
		// floating-point representation
		var cl_LF b = x; // b := x
		var cl_LF sum = cl_float(0,x); // sum := (float 0 x)
		loop {
			var cl_LF new_sum = sum + b;
			if (new_sum == sum) // = sum ?
				break; // ja -> Potenzreihe abbrechen
			sum = new_sum;
			b = (b*a)/(cl_I)((i+1)*(i+2));
			i = i+2;
		}
		powser_value = sum;
	} else {
		// naive3:
		// floating-point representation with smooth precision reduction
		var cl_LF b = x; // b := x
		var cl_LF eps = scale_float(b,-(sintC)d-10);
		var cl_LF sum = cl_float(0,x); // sum := (float 0 x)
		loop {
			var cl_LF new_sum = sum + LF_to_LF(b,actuallen);
			if (new_sum == sum) // = sum ?
				break; // ja -> Potenzreihe abbrechen
			sum = new_sum;
			b = cl_LF_shortenwith(b,eps);
			b = (b*a)/(cl_I)((i+1)*(i+2));
			i = i+2;
		}
		powser_value = sum;
	}
	var cl_LF z = square(powser_value); // sin^2 als Ergebnis
	while (e > e_limit) {
		z = cl_float(1,x) - square(cl_float(1,x) - scale_float(z,1)); // z := 1-(1-2*z)^2
		e--;
	}
	return z;
}}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

// Timings of the three variants, on an i486 33 MHz, running Linux,
// applied to x = sqrt(2)-1 = 0.414...
//   N     naive1  naive2  naive3  ratseries
//    4     0.0064  0.0048  0.0049  0.023
//    6     0.0081  0.0064  0.0065  0.031
//    8     0.0103  0.0085  0.0083  0.038
//   10     0.012   0.011   0.010   0.048
//   25     0.043   0.047   0.035   0.119
//   50     0.15    0.17    0.12    0.37
//  100     0.54    0.67    0.44    1.09
//  250     3.5     4.4     2.8     5.5
//  500    14.7    18.5    11.6    19.4
// 1000    61      78      48      64
// 2500   315     361     243     261
// 2700                   265     270
// 3000                   294     282
// 3500                   339     303
// ==> ratseries faster for N >= 2750.

}  // namespace cln
