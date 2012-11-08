// eulerconst().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/transcendental/cl_LF_tran.h"
#include "float/lfloat/cl_LF.h"
#include "cln/integer.h"
#include "cln/real.h"
#include "base/cl_alloca.h"

namespace cln {

#if 0 // works, but besselintegral4 is always faster

const cl_LF compute_eulerconst_expintegral (uintC len)
{
	// [Jonathan M. Borwein, Peter B. Borwein: Pi and the AGM.
	//  Wiley 1987. Section 10.2.3, exercise 11, p. 336]
	// Use the following formula for the modified exponential integral
	// (valid for Re(z) > 0)
	//   E1(z) := integral(t = z..+infty, exp(-t)/t dt)
	//   E1(z) = - log z - C + sum(n=1..infty, (-1)^(n-1) z^n / (n*n!))
	// [Hint for proving this formula:
	//  1. Learn about the elementary properties of the Gamma function.
	//  2. -C = derivative of Gamma at 1
	//        = lim_{z -> 0} (Gamma(z) - 1/z)
	//        = integral(t = 0..1, (exp(-t)-1)/t dt)
	//          + integral(t = 1..infty, exp(-t)/t dt)
	//  3. Add
	//     0 = integral(t=0..1, 1/(t+1)) - integral(t=1..infty, 1/t(t+1) dt)
	//     to get
	//     -C = integral(t = 0..infty, (exp(-t)/t - 1/t(t+1)) dt)
	//  4. Compute E1(z) + C and note that E1(z) + C + log z is the integral
	//     of an entire function, hence an entire function as well.]
	// Of course we also have the estimate
	//   |E1(z)| < exp(-Re(z)).
	// This means that we can get C by computing
	//   sum(n=1..infty, (-1)^(n-1) z^n / (n*n!)) - log z
	// for large z.
	// In order to get M bits of precision, we first choose z (real) such
	// that exp(-z) < 2^-M. This will make |E1(z)| small enough. z should
	// be chosen as an integer, this is the key to computing the series
	// sum very fast. z = M*log(2) + O(1).
	// Then we choose the number N of terms:
	//   Note than the n-th term's absolute value is (logarithmically)
	//     n*log(z) - n*log(n) + n - 3/2*log(n) - log(sqrt(2 pi)) + o(1).
	//   The derivative of this with respect to n is
	//     log(z) - log(n) - 3/(2n) + o(1/n),
	//   hence is increasing for n < z and decreasing for n > z.
	//   The maximum value is attained at n = z + O(1), and is z + O(log z),
	//   which means that we need z/log(2) + O(log z) bits before the
	//   decimal point.
	//   We can cut off the series when
	//    n*log(z) - n*log(n) + n - 3/2*log(n) - log(sqrt(2 pi)) < -M*log(2)
	//   This happens at n = alpha*z - 3/(2*log(alpha))*log(z) + O(1),
	//   where alpha = 3.591121477... is the solution of
	//     -alpha*log(alpha) + alpha + 1 = 0.
	//   [Use the Newton iteration  alpha --> (alpha+1)/log(alpha)  to
	//    compute this number.]
	// Finally we compute the series's sum as
	//   sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with a(n) = 1, b(n) = n+1, p(n) = z for n=0, -z for n>0, q(n) = n+1.
	// If we computed this with floating-point numbers, we would have
	// to more than double the floating-point precision because of the large
	// extinction which takes place. But luckily we compute with integers.
	var uintC actuallen = len+1; // 1 guard digit
	var uintC z = (uintC)(0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(3.591121477*z);
	CL_ALLOCA_STACK;
	var cl_I* bv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* pv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	for (n = 0; n < N; n++) {
		init1(cl_I, bv[n]) (n+1);
		init1(cl_I, pv[n]) (n==0 ? (cl_I)z : -(cl_I)z);
		init1(cl_I, qv[n]) (n+1);
	}
	var cl_pqb_series series;
	series.bv = bv;
	series.pv = pv; series.qv = qv;
	var cl_LF fsum = eval_rational_series<false>(N,series,actuallen);
	for (n = 0; n < N; n++) {
		bv[n].~cl_I();
		pv[n].~cl_I();
		qv[n].~cl_I();
	}
	fsum = fsum - ln(cl_I_to_LF(z,actuallen)); // log(z) subtrahieren
	return shorten(fsum,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

#endif

#if 0 // works, but besselintegral1 is twice as fast

const cl_LF compute_eulerconst_expintegral1 (uintC len)
{
	// Define f(z) := sum(n=0..infty, z^n/n!) = exp(z)
	// and    g(z) := sum(n=0..infty, H_n*z^n/n!)
	// where H_n := 1/1 + 1/2 + ... + 1/n.
	// The following formula can be proved:
	//   g'(z) - g(z) = (exp(z)-1)/z,
	//   g(z) = exp(z)*(log(z) + c3 + integral(t=..z, exp(-t)/t dt))
	// The Laplace method for determining the asymptotics of an integral
	// or sum yields for real x>0 (the terms n = x+O(x^(1/2)) are
	// dominating):
	//   f(x) = exp(x)*(1 + O(x^(-1/2)))
	//   g(x) = exp(x)*(log(x) + C + O(log(x)*x^(-1/2)))
	// Hence
	//   g(x)/f(x) - log(x) - C = O(log(x)*x^(-1/2))
	// This determines the constant c3, we thus have
	//   g(z) = exp(z)*(log(z) + C + integral(t=z..infty, exp(-t)/t dt))
	// Hence we have for x -> infty:
	//   g(x)/f(x) - log(x) - C == O(exp(-x))
	// This means that we can get C by computing
	//   g(x)/f(x) - log(x)
	// for large x.
	// In order to get M bits of precision, we first choose x (real) such
	// that exp(-x) < 2^-M. This will make the absolute value of the
	// integral small enough. x should be chosen as an integer, this is
	// the key to computing the series sum very fast. x = M*log(2) + O(1).
	// Then we choose the number N of terms:
	//   Note than the n-th term's absolute value is (logarithmically)
	//     n*log(x) - n*log(n) + n - 1/2*log(n) - 1/2*log(2 pi) + o(1).
	//   The derivative of this with respect to n is
	//     log(x) - log(n) - 1/2n + o(1/n),
	//   hence is increasing for n < x and decreasing for n > x.
	//   The maximum value is attained at n = x + O(1), and is
	//   x + O(log x), which means that we need x/log(2) + O(log x)
	//   bits before the decimal point. This also follows from the
	//   asymptotic estimate for f(x).
	//   We can cut off the series when the relative error is < 2^-M,
	//   i.e. when the absolute error is < 2^-M*exp(x), i.e.
	//     n*log(x) - n*log(n) + n - 1/2*log(n) - 1/2*log(2 pi) <
	//     < -M*log(2) + x
	//   This happens at n = e*x - 1/2*log(x) + O(1).
	// Finally we compute the sums of the series f(x) and g(x) with N terms
	// each.
	// We compute f(x) classically and g(x) using the partial sums of f(x).
	var uintC actuallen = len+2; // 2 guard digits
	var uintC x = (uintC)(0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(2.718281828*x);
	var cl_LF one = cl_I_to_LF(1,actuallen);
	var cl_LF fterm = one;
	var cl_LF fsum = fterm;
	var cl_LF gterm = cl_I_to_LF(0,actuallen);
	var cl_LF gsum = gterm;
	var uintC n;
	// After n loops
	//   fterm = x^n/n!, fsum = 1 + x/1! + ... + x^n/n!,
	//   gterm = H_n*x^n/n!, gsum = H_1*x/1! + ... + H_n*x^n/n!.
	for (n = 1; n < N; n++) {
		fterm = The(cl_LF)(fterm*x)/n;
		gterm = (The(cl_LF)(gterm*x) + fterm)/n;
		if (len < 10 || n <= x) {
			fsum = fsum + fterm;
			gsum = gsum + gterm;
		} else {
			// For n > x, the terms are decreasing.
			// So we can reduce the precision accordingly.
			fterm = cl_LF_shortenwith(fterm,one);
			gterm = cl_LF_shortenwith(gterm,one);
			fsum = fsum + LF_to_LF(fterm,actuallen);
			gsum = gsum + LF_to_LF(gterm,actuallen);
		}
	}
	var cl_LF result = gsum/fsum - ln(cl_I_to_LF(x,actuallen));
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(N^2).

#endif

#if 0 // works, but besselintegral4 is always faster

// Same algorithm as expintegral1, but using binary splitting to evaluate
// the sums.
const cl_LF compute_eulerconst_expintegral2 (uintC len)
{
	var uintC actuallen = len+2; // 2 guard digits
	var uintC x = (uintC)(0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(2.718281828*x);
	CL_ALLOCA_STACK;
	var cl_pqd_series_term* args = (cl_pqd_series_term*) cl_alloca(N*sizeof(cl_pqd_series_term));
	var uintC n;
	for (n = 0; n < N; n++) {
		init1(cl_I, args[n].p) (x);
		init1(cl_I, args[n].q) (n+1);
		init1(cl_I, args[n].d) (n+1);
	}
	var cl_pqd_series_result sums;
	eval_pqd_series_aux(N,args,sums);
	// Instead of computing  fsum = 1 + T/Q  and  gsum = V/(D*Q)
	// and then dividing them, to compute  gsum/fsum, we save two
	// divisions by computing  V/(D*(Q+T)).
	var cl_LF result =
	  cl_I_to_LF(sums.V,actuallen)
	  / The(cl_LF)(sums.D * cl_I_to_LF(sums.Q+sums.T,actuallen))
	  - ln(cl_I_to_LF(x,actuallen));
	for (n = 0; n < N; n++) {
		args[n].p.~cl_I();
		args[n].q.~cl_I();
		args[n].d.~cl_I();
	}
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

#endif

// cl_LF compute_eulerconst_besselintegral (uintC len)
	// This is basically the algorithm used in Pari.
	// Define f(z) := sum(n=0..infty, z^n/n!^2)
	// and    g(z) := sum(n=0..infty, H_n*z^n/n!^2)
	// where H_n := 1/1 + 1/2 + ... + 1/n.
	// [f(z) and g(z) are intimately related to the Bessel functions:
	//  f(x^2) = I_0(2*x), g(x^2) = K_0(2*x) + I_0(2*x*)*(C + log(x)).]
	// The following formulas can be proved:
	//   z f''(z) + f'(z) - f(z) = 0,
	//   z g''(z) + g'(z) - g(z) = f'(z),
	//   g(z) = (log(z)/2 + c3 - 1/2 integral(t=..z, 1/(t f(t)^2) dt)) f(z)
	// The Laplace method for determining the asymptotics of an integral
	// or sum yields for real x>0 (the terms n = sqrt(x)+O(x^(1/4)) are
	// dominating):
	//   f(x) = exp(2*sqrt(x))*x^(-1/4)*1/(2*sqrt(pi))*(1 + O(x^(-1/4)))
	//   g(x) = exp(2*sqrt(x))*x^(-1/4)*1/(2*sqrt(pi))*
	//          (1/2*log(x) + C + O(log(x)*x^(-1/4)))
	// Hence
	//   g(x)/f(x) - 1/2*log(x) - C = O(log(x)*x^(-1/4))
	// This determines the constant c3, we thus have
	// g(z)= (log(z)/2 + C + 1/2 integral(t=z..infty, 1/(t f(t)^2) dt)) f(z)
	// Hence we have for x -> infty:
	//   g(x)/f(x) - 1/2*log(x) - C == pi*exp(-4*sqrt(x)) [approx.]
	// This means that we can get C by computing
	//   g(x)/f(x) - 1/2*log(x)
	// for large x.
	// In order to get M bits of precision, we first choose x (real) such
	// that exp(-4*sqrt(x)) < 2^-(M+2). This will make the absolute value
	// of the integral small enough. x should be chosen as an integer,
	// this is the key to computing the series sum very fast. sqrt(x)
	// need not be an integer. Set sx = sqrt(x).
	// sx = M*log(2)/4 + O(1).
	// Then we choose the number N of terms:
	//   Note than the n-th term's absolute value is (logarithmically)
	//     n*log(x) - 2*n*log(n) + 2*n - log(n) - log(2 pi) + o(1).
	//   The derivative of this with respect to n is
	//     log(x) - 2*log(n) - 1/n + o(1/n),
	//   hence is increasing for n < sx and decreasing for n > sx.
	//   The maximum value is attained at n = sx + O(1), and is
	//   2*sx + O(log x), which means that we need 2*sx/log(2) + O(log x)
	//   bits before the decimal point. This also follows from the
	//   asymptotic estimate for f(x).
	//   We can cut off the series when the relative error is < 2^-M,
	//   i.e. when the absolute error is
	//      < 2^-M*exp(2*sx)*sx^(-1/2)*1/(2*sqrt(pi)),
	//   i.e.
	//     n*log(x) - 2*n*log(n) + 2*n - log(n) - log(2 pi) <
	//     < -M*log(2) + 2*sx - 1/2*log(sx) - log(2 sqrt(pi))
	//   This happens at n = alpha*sx - 1/(4*log(alpha))*log(sx) + O(1),
	//   where alpha = 3.591121477... is the solution of
	//     -alpha*log(alpha) + alpha + 1 = 0.
	//   [Use the Newton iteration  alpha --> (alpha+1)/log(alpha)  to
	//    compute this number.]
	// Finally we compute the sums of the series f(x) and g(x) with N terms
	// each.

const cl_LF compute_eulerconst_besselintegral1 (uintC len)
{
	// We compute f(x) classically and g(x) using the partial sums of f(x).
	var uintC actuallen = len+1; // 1 guard digit
	var uintC sx = (uintC)(0.25*0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(3.591121477*sx);
	var cl_I x = square((cl_I)sx);
	var cl_LF eps = scale_float(cl_I_to_LF(1,LF_minlen),-(sintC)(sx*2.88539+10));
	var cl_LF fterm = cl_I_to_LF(1,actuallen);
	var cl_LF fsum = fterm;
	var cl_LF gterm = cl_I_to_LF(0,actuallen);
	var cl_LF gsum = gterm;
	var uintC n;
	// After n loops
	//   fterm = x^n/n!^2, fsum = 1 + x/1!^2 + ... + x^n/n!^2,
	//   gterm = H_n*x^n/n!^2, gsum = H_1*x/1!^2 + ... + H_n*x^n/n!^2.
	for (n = 1; n < N; n++) {
		fterm = The(cl_LF)(fterm*x)/square((cl_I)n);
		gterm = (The(cl_LF)(gterm*x)/(cl_I)n + fterm)/(cl_I)n;
		if (len < 10 || n <= sx) {
			fsum = fsum + fterm;
			gsum = gsum + gterm;
		} else {
			// For n > sx, the terms are decreasing.
			// So we can reduce the precision accordingly.
			fterm = cl_LF_shortenwith(fterm,eps);
			gterm = cl_LF_shortenwith(gterm,eps);
			fsum = fsum + LF_to_LF(fterm,actuallen);
			gsum = gsum + LF_to_LF(gterm,actuallen);
		}
	}
	var cl_LF result = gsum/fsum - ln(cl_I_to_LF(sx,actuallen));
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(N^2).

#if 0 // works, but besselintegral1 is faster

const cl_LF compute_eulerconst_besselintegral2 (uintC len)
{
	// We compute the sum of the series f(x) as
	//   sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with a(n) = 1, b(n) = 1, p(n) = x for n>0, q(n) = n^2 for n>0.
	// and the sum of the series g(x) as
	//   sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with
	// a(n) = HN_{n+1}, b(n) = 1, p(n) = x, q(n) = (n+1)^2 * HD_{n+1}/HD_{n}
	// where HD_n := lcm(1,...,n) and HN_n := HD_n * H_n. (Note that
	// HD_n need not be the lowest possible denominator of H_n. For
	// example, n=6: H_6 = 49/20, but HD_6 = 60.)
	// WARNING: The memory used by this algorithm grown quadratically in N.
	// (Because HD_n grows like exp(n), hence HN_n grows like exp(n) as
	// well, and we store all HN_n values in an array!)
	var uintC actuallen = len+1; // 1 guard digit
	var uintC sx = (uintC)(0.25*0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(3.591121477*sx);
	var cl_I x = square((cl_I)sx);
	CL_ALLOCA_STACK;
	var cl_I* av = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* pv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	// Evaluate f(x).
	init1(cl_I, pv[0]) (1);
	init1(cl_I, qv[0]) (1);
	for (n = 1; n < N; n++) {
		init1(cl_I, pv[n]) (x);
		init1(cl_I, qv[n]) ((cl_I)n*(cl_I)n);
	}
	var cl_pq_series fseries;
	fseries.pv = pv; fseries.qv = qv;
	var cl_LF fsum = eval_rational_series<false>(N,fseries,actuallen);
	for (n = 0; n < N; n++) {
		pv[n].~cl_I();
		qv[n].~cl_I();
	}
	// Evaluate g(x).
	var cl_I HN = 0;
	var cl_I HD = 1;
	for (n = 0; n < N; n++) {
		// Now HN/HD = H_n.
		var cl_I Hu = gcd(HD,n+1);
		var cl_I Hv = exquopos(n+1,Hu);
		HN = HN*Hv + exquopos(HD,Hu);
		HD = HD*Hv;
		// Now HN/HD = H_{n+1}.
		init1(cl_I, av[n]) (HN);
		init1(cl_I, pv[n]) (x);
		init1(cl_I, qv[n]) (Hv*(cl_I)(n+1)*(cl_I)(n+1));
	}
	var cl_pqa_series gseries;
	gseries.av = av;
	gseries.pv = pv; gseries.qv = qv;
	var cl_LF gsum = eval_rational_series<false>(N,gseries,actuallen);
	for (n = 0; n < N; n++) {
		av[n].~cl_I();
		pv[n].~cl_I();
		qv[n].~cl_I();
	}
	var cl_LF result = gsum/fsum - ln(cl_I_to_LF(sx,actuallen));
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(N^2).
// Memory consumption: O(N^2).

// Same algorithm as besselintegral2, but without quadratic memory consumption.
#define cl_rational_series_for_g cl_rational_series_for_besselintegral3_g
struct cl_rational_series_for_g : cl_pqa_series_stream {
	uintL n;
	cl_I HN;
	cl_I HD;
	cl_I x;
	static cl_pqa_series_term computenext (cl_pqa_series_stream& thisss)
	{
		var cl_rational_series_for_g& thiss = (cl_rational_series_for_g&)thisss;
		var uintL n = thiss.n;
		// Now HN/HD = H_n.
		var cl_I Hu = gcd(thiss.HD,n+1);
		var cl_I Hv = exquopos(n+1,Hu);
		thiss.HN = thiss.HN*Hv + exquopos(thiss.HD,Hu);
		thiss.HD = thiss.HD*Hv;
		// Now HN/HD = H_{n+1}.
		var cl_pqa_series_term result;
		result.p = thiss.x;
		result.q = Hv*(cl_I)(n+1)*(cl_I)(n+1);
		result.a = thiss.HN;
		thiss.n = n+1;
		return result;
	}
	cl_rational_series_for_g (const cl_I& _x)
		: cl_pqa_series_stream (cl_rational_series_for_g::computenext),
		  n (0), HN (0), HD (1), x (_x) {}
};
const cl_LF compute_eulerconst_besselintegral3 (uintC len)
{
	var uintC actuallen = len+1; // 1 guard digit
	var uintC sx = (uintC)(0.25*0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(3.591121477*sx);
	var cl_I x = square((cl_I)sx);
	CL_ALLOCA_STACK;
	var cl_I* pv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	// Evaluate f(x).
	init1(cl_I, pv[0]) (1);
	init1(cl_I, qv[0]) (1);
	for (n = 1; n < N; n++) {
		init1(cl_I, pv[n]) (x);
		init1(cl_I, qv[n]) ((cl_I)n*(cl_I)n);
	}
	var cl_pq_series fseries;
	fseries.pv = pv; fseries.qv = qv;
	var cl_LF fsum = eval_rational_series<false>(N,fseries,actuallen);
	for (n = 0; n < N; n++) {
		pv[n].~cl_I();
		qv[n].~cl_I();
	}
	// Evaluate g(x).
	var cl_rational_series_for_g gseries = cl_rational_series_for_g(x);
	var cl_LF gsum = eval_rational_series<false>(N,gseries,actuallen);
	var cl_LF result = gsum/fsum - ln(cl_I_to_LF(sx,actuallen));
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(N^2).

#endif

// Same algorithm as besselintegral1, but using binary splitting to evaluate
// the sums.
const cl_LF compute_eulerconst_besselintegral4 (uintC len)
{
	var uintC actuallen = len+2; // 2 guard digits
	var uintC sx = (uintC)(0.25*0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(3.591121477*sx);
	var cl_I x = square(cl_I(sx));
	struct rational_series_stream : cl_pqd_series_stream {
		uintC n;
		cl_I x;
		static cl_pqd_series_term computenext (cl_pqd_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_pqd_series_term result;
			result.p = thiss.x;
			result.q = square(cl_I(n+1));
			result.d = n+1;
			thiss.n = n+1;
			return result;
		}
		rational_series_stream (uintC n_, const cl_I& x_)
			: cl_pqd_series_stream (rational_series_stream::computenext),
			  n (n_), x (x_) {}
	} series(0,x);
	var cl_pqd_series_result<cl_R> sums;
	eval_pqd_series_aux(N,series,sums,actuallen);
	// Instead of computing  fsum = 1 + T/Q  and  gsum = V/(D*Q)
	// and then dividing them, to compute  gsum/fsum, we save two
	// divisions by computing  V/(D*(Q+T)).
	var cl_LF result =
	  cl_R_to_LF(sums.V,actuallen)
	  / The(cl_LF)(sums.D * cl_R_to_LF(sums.Q+sums.T,actuallen))
	  - ln(cl_R_to_LF(sx,actuallen));
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an i486 33 MHz, running Linux.
//    N      exp      exp1    exp2    bessel1   bessel2   bessel3   bessel4
//    10     0.51     0.28    0.52     0.11      0.16      0.16      0.15
//    25     2.23     0.83    2.12     0.36      0.62      0.63      0.62
//    50     6.74     2.23    6.54     0.95      1.95      1.97      1.95
//   100    19.1      6.74   20.6      2.96      6.47      6.42      6.3
//   250    84       37.4    78       16.3      33.6      32.0      28.8
//   500   230      136.5   206       60.5       ---     111        85
//  1000   591      520     536      229         ---     377       241
//  1050                             254                           252
//  1100                             277                           266
//  2500  1744             2108    (1268)                          855 (run)
//  2500  1845             2192    (1269)                          891 (real)
//
// asymp.  FAST      N^2     FAST    N^2        N^2      N^2       FAST
// (FAST means O(log(N)^2*M(N)))
//
// The break-even point between "bessel1" and "bessel4" is at about N = 1050.
const cl_LF compute_eulerconst (uintC len)
{
	if (len >= 1050)
		return compute_eulerconst_besselintegral4(len);
	else
		return compute_eulerconst_besselintegral1(len);
}

const cl_LF eulerconst (uintC len)
{
	var uintC oldlen = TheLfloat(cl_LF_eulerconst())->len; // vorhandene Länge
	if (len < oldlen)
		return shorten(cl_LF_eulerconst(),len);
	if (len == oldlen)
		return cl_LF_eulerconst();

	// TheLfloat(cl_LF_eulerconst())->len um mindestens einen konstanten Faktor
	// > 1 wachsen lassen, damit es nicht zu häufig nachberechnet wird:
	var uintC newlen = len;
	oldlen += floor(oldlen,2); // oldlen * 3/2
	if (newlen < oldlen)
		newlen = oldlen;

	// gewünschte > vorhandene Länge -> muß nachberechnen:
	cl_LF_eulerconst() = compute_eulerconst(newlen);
	return (len < newlen ? shorten(cl_LF_eulerconst(),len) : cl_LF_eulerconst());
}

}  // namespace cln
