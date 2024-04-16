// pi().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/transcendental/cl_LF_tran.h"
#include "float/lfloat/cl_LF.h"
#include "cln/integer.h"
#include "base/cl_alloca.h"

namespace cln {

ALL_cl_LF_OPERATIONS_SAME_PRECISION()

// For the next algorithms, I warmly recommend
// [Jörg Arndt: hfloat documentation, august 1996,
//  http://www.tu-chemnitz.de/~arndt/ hfdoc.dvi
//  But beware of the typos in his formulas! ]

const cl_LF compute_pi_brent_salamin (uintC len)
{
	// Methode:
	// [Richard P. Brent: Fast multiple-precision evaluation of elementary
	//  functions. J. ACM 23(1976), 242-251.]
	// [Jonathan M. Borwein, Peter B. Borwein: Pi and the AGM.
	//  Wiley 1987. Algorithm 2.2, p. 48.]
	// [Jörg Arndt, formula (4.51)-(4.52).]
	//    pi = AGM(1,1/sqrt(2))^2 * 2/(1 - sum(k=0..infty, 2^k c_k^2)).
	// where the AGM iteration reads
	//    a_0 := 1, b_0 := 1/sqrt(2).
	//    a_(k+1) := (a_k + b_k)/2, b_(k+1) := sqrt(a_k*b_k).
	// and we set
	//    c_k^2 := a_k^2 - b_k^2,
	// i.e. c_0^2 = 1/2,
	//      c_(k+1)^2 = a_(k+1)^2 - b_(k+1)^2 = (a_k - b_k)^2/4
	//                = (a_k - a_(k+1))^2.
	// (Actually c_(k+1) is _defined_ as  = a_k - a_(k+1) = (a_k - b_k)/2.)
	// d=len, n:=intDsize*d. Verwende Long-Floats mit intDsize*(d+1)
	// Mantissenbits.
	// (let* ((a (coerce 1 'long-float)) ; 1
	//        (b (sqrt (scale-float a -1))) ; 2^-(1/2)
	//        (eps (scale-float a (- n))) ; 2^-n
	//        (t (scale-float a -2)) ; 1/4
	//        (k 0)
	//       )
	//   (loop
	//     (when (< (- a b) eps) (return))
	//     (let ((y a))
	//       (setq a (scale-float (+ a b) -1))
	//       (setq b (sqrt (* b y)))
	//       (setq t (- t (scale-float (expt (- a y) 2) k)))
	//     )
	//     (incf k)
	//   )
	//   (/ (expt a 2) t)
	// )
	var uintC actuallen = len + 1; // 1 guard digit
	var uintE uexp_limit = LF_exp_mid - intDsize*len;
	// Ein Long-Float ist genau dann betragsmäßig <2^-n, wenn
	// sein Exponent < LF_exp_mid-n = uexp_limit ist.
	var cl_LF a = cl_I_to_LF(1,actuallen);
	var cl_LF b = sqrt(scale_float(a,-1));
	var uintL k = 0;
	var cl_LF t = scale_float(a,-2);
	until (TheLfloat(a-b)->expo < uexp_limit) {
		// |a-b| < 2^-n -> fertig
		var cl_LF new_a = scale_float(a+b,-1); // (a+b)/2
		b = sqrt(a*b);
		var cl_LF a_diff = new_a - a;
		t = t - scale_float(square(a_diff),k);
		a = new_a;
		k++;
	}
	var cl_LF pires = square(a)/t; // a^2/t
	return shorten(pires,len); // verkürzen und fertig
}
// Bit complexity (N := len): O(log(N)*M(N)).

const cl_LF compute_pi_brent_salamin_quartic (uintC len)
{
	// See [Borwein, Borwein, section 1.4, exercise 3, p. 17].
	// See also [Jörg Arndt], formulas (4.52) and (3.30)[wrong!].
	// As above, we are using the formula
	//    pi = AGM(1,1/sqrt(2))^2 * 2/(1 - sum(k=0..infty, 2^k c_k^2)).
	// where the AGM iteration reads
	//    a_0 := 1, b_0 := 1/sqrt(2).
	//    a_(k+1) := (a_k + b_k)/2, b_(k+1) := sqrt(a_k*b_k).
	// But we keep computing with
	//    wa_k := sqrt(a_k)  and  wb_k := sqrt(b_k)
	// at the same time and do two iteration steps at once.
	// The iteration takes now takes the form
	//    wa_0 := 1, wb_0 := 2^-1/4,
	//    (wa_k, wb_k, a_k, b_k)
	//    -->         ((wa_k^2 + wb_k^2)/2, wa_k*wb_k)
	//    -->         (((wa_k + wb_k)/2)^2, sqrt(wa_k*wb_k*(wa_k^2 + wb_k^2)/2)),
	//    i.e. wa_(k+2) = (wa_k + wb_k)/2 and
	//         wb_(k+2) = sqrt(sqrt(wa_k*wb_k*(wa_k^2 + wb_k^2)/2)).
	// For the sum, we can group two successive items together:
	//   2^k * c_k^2 + 2^(k+1) * c_(k+1)^2
	//   = 2^k * [(a_k^2 - b_k^2) + 2*((a_k - b_k)/2)^2]
	//   = 2^k * [3/2*a_k^2 - a_k*b_k - 1/2*b_k^2]
	//   = 2^k * [2*a_k^2 - 1/2*(a_k+b_k)^2]
	//   = 2^(k+1) * [a_k^2 - ((a_k+b_k)/2)^2]
	//   = 2^(k+1) * [wa_k^4 - ((wa_k^2+wb_k^2)/2)^2].
	// Hence,
	//   pi = AGM(1,1/sqrt(2))^2 * 1/(1/2 - sum(k even, 2^k*[....])).
	var uintC actuallen = len + 1; // 1 guard digit
	var uintE uexp_limit = LF_exp_mid - intDsize*len;
	var cl_LF one = cl_I_to_LF(1,actuallen);
	var cl_LF a = one;
	var cl_LF wa = one;
	var cl_LF b = sqrt(scale_float(one,-1));
	var cl_LF wb = sqrt(b);
	// We keep a = wa^2, b = wb^2.
	var uintL k = 0;
	var cl_LF t = scale_float(one,-1);
	until (TheLfloat(wa-wb)->expo < uexp_limit) {
		// |wa-wb| < 2^-n -> fertig
		var cl_LF wawb = wa*wb;
		var cl_LF new_wa = scale_float(wa+wb,-1);
		var cl_LF a_b = scale_float(a+b,-1);
		var cl_LF new_a = scale_float(a_b+wawb,-1);
		var cl_LF new_b = sqrt(wawb*a_b);
		var cl_LF new_wb = sqrt(new_b);
		t = t - scale_float((a - a_b)*(a + a_b),k);
		a = new_a; wa = new_wa;
		b = new_b; wb = new_wb;
		k += 2;
	}
	var cl_LF pires = square(a)/t;
	return shorten(pires,len); // verkürzen und fertig
}
// Bit complexity (N := len): O(log(N)*M(N)).

const cl_LF compute_pi_ramanujan_163 (uintC len)
{
	// 1/pi = 1/sqrt(-1728 J)
	//        * sum(n=0..infty, (6n)! (A+nB) / 12^(3n) (3n)! n!^3 J^n)
	// mit J = -53360^3 = - (2^4 5 23 29)^3
	//     A = 163096908 = 2^2 3 13 1045493
	//     B = 6541681608 = 2^3 3^3 7 11 19 127 163
	// See [Jörg Arndt], formulas (4.27)-(4.30).
	// This is also the formula used in Pari.
	// The absolute value of the n-th summand is approximately
	//   |J|^-n * n^(-1/2) * B/(2*pi^(3/2)),
	// hence every summand gives more than 14 new decimal digits
	// in precision.
	// The sum is best evaluated using fixed-point arithmetic,
	// so that the precision is reduced for the later summands.
	var uintC actuallen = len + 4; // 4 guard digits
	var sintC scale = intDsize*actuallen;
	static const cl_I A = "163096908";
	static const cl_I B = "6541681608";
	//static const cl_I J1 = "10939058860032000"; // 72*abs(J)
	static const cl_I J2 = "333833583375"; // odd part of J1
	var cl_I sum = 0;
	var cl_I n = 0;
	var cl_I factor = ash(1,scale);
	while (!zerop(factor)) {
		sum = sum + factor * (A+n*B);
		factor = factor * ((6*n+1)*(2*n+1)*(6*n+5));
		n = n+1;
		factor = truncate1(factor,n*n*n*J2);
		// Finally divide by 2^15 and change sign.
		if (minusp(factor))
			factor = (-factor) >> 15;
		else
			factor = -(factor >> 15);
	}
	var cl_LF fsum = scale_float(cl_I_to_LF(sum,actuallen),-scale);
	static const cl_I J3 = "262537412640768000"; // -1728*J
	var cl_LF pires = sqrt(cl_I_to_LF(J3,actuallen)) / fsum;
	return shorten(pires,len); // verkürzen und fertig
}
// Bit complexity (N := len): O(N^2).

const cl_LF compute_pi_ramanujan_163_fast (uintC len)
{
	// Same formula as above, using a binary splitting evaluation.
	// See [Borwein, Borwein, section 10.2.3].
	struct rational_series_stream : cl_pqa_series_stream {
		uintC n;
		static cl_pqa_series_term computenext (cl_pqa_series_stream& thisss)
		{
			static const cl_I A = "163096908";
			static const cl_I B = "6541681608";
			static const cl_I J1 = "10939058860032000"; // 72*abs(J)
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_pqa_series_term result;
			if (n==0) {
				result.p = 1;
				result.q = 1;
			} else {
				result.p = -((cl_I)(6*n-5)*(cl_I)(2*n-1)*(cl_I)(6*n-1));
				result.q = (cl_I)n*(cl_I)n*(cl_I)n*J1;
			}
			result.a = A+n*B;
			thiss.n = n+1;
			return result;
		}
		rational_series_stream ()
			: cl_pqa_series_stream (rational_series_stream::computenext),
			  n (0) {}
	} series;
	var uintC actuallen = len + 4; // 4 guard digits
	static const cl_I A = "163096908";
	static const cl_I B = "6541681608";
	static const cl_I J1 = "10939058860032000"; // 72*abs(J)
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = A+n*B, b(n) = 1,
	//   p(n) = -(6n-5)(2n-1)(6n-1) for n>0,
	//   q(n) = 72*|J|*n^3 for n>0.
	var const uintC n_slope = (uintC)(intDsize*32*0.02122673)+1;
	// n_slope >= 32*intDsize*log(2)/log(|J|), normally n_slope = 22.
	var uintC N = (n_slope*actuallen)/32 + 1;
	// N > intDsize*log(2)/log(|J|) * actuallen, hence
	// |J|^-N < 2^(-intDsize*actuallen).
	var cl_LF fsum = eval_rational_series<false>(N,series,actuallen,actuallen);
	static const cl_I J3 = "262537412640768000"; // -1728*J
	var cl_LF pires = sqrt(cl_I_to_LF(J3,actuallen)) / fsum;
	return shorten(pires,len); // verkürzen und fertig
}
// Bit complexity (N := len): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an i486 33 MHz, running Linux.
//    N      Brent   Brent4  R 163   R 163 fast
//    10     0.0079  0.0079  0.0052  0.0042
//    25     0.026   0.026   0.014   0.012
//    50     0.085   0.090   0.037   0.033
//   100     0.29    0.29    0.113   0.098
//   250     1.55    1.63    0.60    0.49
//   500     5.7     5.7     2.24    1.71
//  1000    21.6    22.9     8.5     5.5
//  2500    89      95      49      19.6
//  5000   217     218     188      49
// 10000   509     540     747     117
// 25000  1304    1310    4912     343
// We see that
//   - "Brent4" isn't worth it: No speed improvement over "Brent".
//   - "R 163" is pretty fast at the beginning, but it is an O(N^2)
//     algorithm, hence it loses in the end,
//   - "R 163 fast", which uses the same formula as "R 163", but evaluates
//     it using binary splitting, is an O(log N * M(N)) algorithm, and
//     outperforms all of the others.

const cl_LF pi (uintC len)
{
	var uintC oldlen = TheLfloat(cl_LF_pi())->len; // vorhandene Länge
	if (len < oldlen)
		return shorten(cl_LF_pi(),len);
	if (len == oldlen)
		return cl_LF_pi();

	// TheLfloat(cl_LF_pi())->len um mindestens einen konstanten Faktor
	// > 1 wachsen lassen, damit es nicht zu häufig nachberechnet wird:
	var uintC newlen = len;
	oldlen += floor(oldlen,2); // oldlen * 3/2
	if (newlen < oldlen)
		newlen = oldlen;

	// gewünschte > vorhandene Länge -> muß nachberechnen:
	cl_LF_pi() = compute_pi_ramanujan_163_fast(newlen);
	return (len < newlen ? shorten(cl_LF_pi(),len) : cl_LF_pi());
}

}  // namespace cln
