// lnx().

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

// cl_F lnx_naive (const cl_F& x)
// cl_LF lnx_naive (const cl_LF& x)
//
// Methode:
// y:=x-1, e := Exponent aus (decode-float y), d := (float-digits y)
// Bei y=0.0 oder e<=-d liefere y
//   (denn bei e<=-d ist y/2 < 2^(-d)/2 = 2^(-d-1), also
//   0 <= y - ln(x) < y^2/2 < 2^(-d-1)*y
//   also ist ln(x)/y, auf d Bits gerundet, gleich y).
// Bei e<=-sqrt(d) verwende die Potenzreihe
//   ln(x) = sum(j=0..inf,(-1)^j*y^(j+1)/(j+1)):
//   a:=-y, b:=y, i:=1, sum:=0,
//   while (/= sum (setq sum (+ sum (/ b i)))) do i:=i+1, b:=b*a.
//   Ergebnis sum.
// Sonst setze y := sqrt(x), berechne rekursiv z:=ln(y)
//   und liefere 2*z = (scale-float z 1).
// Aufwand: asymptotisch d^0.5*M(d) = d^2.5 .

const cl_LF lnx_naive (const cl_LF& x)
{
	var cl_LF y = x-cl_float(1,x);
	if (zerop_inline(y)) // y=0.0 -> y als Ergebnis
		return y;
	var uintC actuallen = TheLfloat(x)->len;
	var uintC d = float_digits(x);
	var sintE e = float_exponent_inline(y);
	if (e <= -(sintC)d) // e <= -d ?
		return y; // ja -> y als Ergebnis
 {	Mutable(cl_LF,x);
	var uintL k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-limit_slope*floor(sqrt(d)) kann die Potenzreihe
	// angewandt werden.
	// Wähle für ln(1+y), naive1: limit_slope = 1.0,
	//       für ln(1+y), naive2: limit_slope = 11/16 = 0.7,
	//       für atanh(z), naive1: limit_slope = 0.6,
	//       für atanh(z), naive1: limit_slope = 0.5.
	var sintL e_limit = -1-floor(isqrtC(d),2); // -1-floor(sqrt(d))
	while (e > e_limit) {
		// e > -1-floor(sqrt(d)) -> muß |y| verkleinern.
		x = sqrt(x); // x := (sqrt x)
		y = x-cl_float(1,x); // y := (- x 1) und
		e = float_exponent_inline(y); // e neu berechnen
		k = k+1; // k:=k+1
	}
	if (0) {
		// Potenzreihe ln(1+y) anwenden:
		var int i = 1;
		var cl_LF sum = cl_float(0,x); // sum := (float 0 x)
		var cl_LF a = -y;
		var cl_LF b = y;
		if (0) {
			// naive1:
			// floating-point representation
			loop {
				var cl_LF new_sum = sum + b/(cl_I)i; // (+ sum (/ b i))
				if (new_sum == sum) // = sum ?
					break; // ja -> Potenzreihe abbrechen
				sum = new_sum;
				b = b*a;
				i = i+1;
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
				i = i+1;
			}
		}
		return scale_float(sum,k); // sum als Ergebnis, wegen Rekursion noch mal 2^k
	} else {
		var cl_LF z = y / (x+cl_float(1,x));
		// Potenzreihe atanh(z) anwenden:
		var int i = 1;
		var cl_LF a = square(z); // a = x^2
		var cl_LF b = cl_float(1,x); // b := (float 1 x)
		var cl_LF sum = cl_float(0,x); // sum := (float 0 x)
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
		return scale_float(sum*z,k+1); // 2*sum*z als Ergebnis, wegen Rekursion noch mal 2^k
	}
}}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

const cl_F lnx_naive (const cl_F& x)
{
	if (longfloatp(x)) {
		DeclareType(cl_LF,x);
		return lnx_naive(x);
	}
	var cl_F y = x-cl_float(1,x);
	if (zerop(y)) // y=0.0 -> y als Ergebnis
		return y;
	var uintC d = float_digits(x);
	var sintE e = float_exponent(y);
	if (e <= -(sintC)d) // e <= -d ?
		return y; // ja -> y als Ergebnis
 {	Mutable(cl_F,x);
	var uintL k = 0; // Rekursionszähler k:=0
	// Bei e <= -1-floor(sqrt(d)) kann die Potenzreihe angewandt werden.
	var sintL e_limit = -1-isqrtC(d); // -1-floor(sqrt(d))
	while (e > e_limit) {
		// e > -1-floor(sqrt(d)) -> muß |y| verkleinern.
		x = sqrt(x); // x := (sqrt x)
		y = x-cl_float(1,x); // y := (- x 1) und
		e = float_exponent(y); // e neu berechnen
		k = k+1; // k:=k+1
	}
	// Potenzreihe anwenden:
	var int i = 1;
	var cl_F sum = cl_float(0,x); // sum := (float 0 x)
	var cl_F a = -y;
	var cl_F b = y;
	loop {
		var cl_F new_sum = sum + b/(cl_I)i; // (+ sum (/ b i))
		if (new_sum == sum) // = sum ?
			break; // ja -> Potenzreihe abbrechen
		sum = new_sum;
		b = b*a;
		i = i+1;
	}
	return scale_float(sum,k); // sum als Ergebnis, wegen Rekursion noch mal 2^k
}}
// Bit complexity (N = length(x)): O(N^(1/2)*M(N)).

const cl_LF lnx_ratseries (const cl_LF& x)
{
	// Method:
	// Based on the same ideas as expx_ratseries.
	// y := 0.
	// Loop
	//   [x*exp(y) is invariant]
	//   x' := x-1. If x' = 0, terminate the loop.
	//   Choose approximation y' of log(x) = log(1+x'):
	//     If |x'| >= 1/2, set y' = 1/2 * sign(x').
	//     If |x'| < 2^-n with n maximal, set
	//       y' = truncate(x'*2^(2n))/2^(2n).
	//   Set y := y + y' and x := x*exp(-y').
	var uintC len = TheLfloat(x)->len;
 {	Mutable(cl_LF,x);
	var cl_LF y = cl_I_to_LF(0,len);
	loop {
		var cl_LF x1 = x + cl_I_to_LF(-1,len);
		var cl_idecoded_float x1_ = integer_decode_float(x1);
		// x1 = (-1)^sign * 2^exponent * mantissa
		if (zerop(x1_.mantissa))
			break;
		var uintC lm = integer_length(x1_.mantissa);
		var uintE me = cl_I_to_UE(- x1_.exponent);
		var cl_I p;
		var uintE lq;
		var bool last_step = false;
		if (lm >= me) { // |x'| >= 1/2 ?
			p = x1_.sign; // 1 or -1
			lq = 1;
		} else {
			var uintE n = me - lm; // |x'| < 2^-n with n maximal
			// Set p to the first n bits of |x'|:
			if (lm > n) {
				p = x1_.mantissa >> (lm - n);
				lq = 2*n;
			} else {
				p = x1_.mantissa;
				lq = lm + n;
			}
			if (minusp(x1_.sign)) { p = -p; }
			// If 2*n >= lm = intDsize*len, then within our
			// precision exp(-y') = 1-y', (because |y'^2| < 2^-lm),
			// and we know a priori that the iteration will stop
			// after the next big multiplication. This saves one
			// big multiplication at the end.
			if (2*n >= lm)
				last_step = true;
		}
		y = y + scale_float(cl_I_to_LF(p,len),-(sintE)lq);
		if (last_step)
			break;
		x = x * cl_exp_aux(-p,lq,len);
	}
	return y;
}}
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an i486 33 MHz, running Linux,
// applied to x = sqrt(sqrt(2)) = 1.189...
//   N      ln(1+y) ln(1+y) atanh z atanh z   exp
//          naive1  naive2  naive1  naive2  ratseries
//   10     0.019   0.016   0.013   0.012   0.036
//   25     0.077   0.056   0.057   0.040   0.087
//   50     0.30    0.21    0.23    0.15    0.21
//  100     1.24    0.81    0.92    0.59    0.61
//  250     8.8     5.8     6.3     4.3     2.77
//  500    43.9    28.8    29.7    21.0     9.8
// 1000   223     149     144     107      30
// ==> ratseries faster for N >= 110. (N = length before extended by the caller.)

}  // namespace cln
