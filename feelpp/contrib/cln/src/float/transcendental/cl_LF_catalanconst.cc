// catalanconst().

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

const cl_LF compute_catalanconst_ramanujan (uintC len)
{
	// [Jonathan M. Borwein, Peter B. Borwein: Pi and the AGM.
	//  Wiley 1987. Section 11.3, exercise 16 g, p. 386]
	// G = 3/8 * sum(n=0..infty, n!^2 / (2n+1)!*(2n+1))
	//     + pi/8 * log(2+sqrt(3)).
	// Every summand gives 0.6 new decimal digits in precision.
	// The sum is best evaluated using fixed-point arithmetic,
	// so that the precision is reduced for the later summands.
	var uintC actuallen = len + 2; // 2 guard digits
	var sintC scale = intDsize*actuallen;
	var cl_I sum = 0;
	var cl_I n = 0;
	var cl_I factor = ash(1,scale);
	while (!zerop(factor)) {
		sum = sum + truncate1(factor,2*n+1);
		n = n+1;
		factor = truncate1(factor*n,2*(2*n+1));
	}
	var cl_LF fsum = scale_float(cl_I_to_LF(sum,actuallen),-scale);
	var cl_LF g =
	  scale_float(The(cl_LF)(3*fsum)
	              + The(cl_LF)(pi(actuallen))
	                * ln(cl_I_to_LF(2,actuallen)+sqrt(cl_I_to_LF(3,actuallen))),
	              -3);
	return shorten(g,len); // verkürzen und fertig
}
// Bit complexity (N := len): O(N^2).

const cl_LF compute_catalanconst_ramanujan_fast (uintC len)
{
	// Same formula as above, using a binary splitting evaluation.
	// See [Borwein, Borwein, section 10.2.3].
	struct rational_series_stream : cl_pqb_series_stream {
		cl_I n;
		static cl_pqb_series_term computenext (cl_pqb_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var cl_I n = thiss.n;
			var cl_pqb_series_term result;
			if (n==0) {
				result.p = 1;
				result.q = 1;
				result.b = 1;
			} else {
				result.p = n;
				result.b = 2*n+1;
				result.q = result.b << 1; // 2*(2*n+1)
			}
			thiss.n = n+1;
			return result;
		}
		rational_series_stream ()
			: cl_pqb_series_stream (rational_series_stream::computenext),
			  n (0) {}
	} series;
	var uintC actuallen = len + 2; // 2 guard digits
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = 1, b(n) = 2*n+1,
	//   p(n) = n for n>0, q(n) = 2*(2*n+1) for n>0.
	var uintC N = (intDsize/2)*actuallen;
	// 4^-N <= 2^(-intDsize*actuallen).
	var cl_LF fsum = eval_rational_series<false>(N,series,actuallen,actuallen);
	var cl_LF g =
	  scale_float(The(cl_LF)(3*fsum)
	              + The(cl_LF)(pi(actuallen))
	                * ln(cl_I_to_LF(2,actuallen)+sqrt(cl_I_to_LF(3,actuallen))),
	              -3);
	return shorten(g,len); // verkürzen und fertig
}
// Bit complexity (N := len): O(log(N)^2*M(N)).

const cl_LF compute_catalanconst_expintegral1 (uintC len)
{
	// We compute f(x) classically and g(x) using the partial sums of f(x).
	var uintC actuallen = len+2; // 2 guard digits
	var uintC x = (uintC)(0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(2.718281828*x);
	var cl_LF fterm = cl_I_to_LF(1,actuallen);
	var cl_LF fsum = fterm;
	var cl_LF gterm = fterm;
	var cl_LF gsum = gterm;
	var uintC n;
	// After n loops
	//   fterm = x^n/n!, fsum = 1 + x/1! + ... + x^n/n!,
	//   gterm = S_n*x^n/n!, gsum = S_0*x^0/0! + ... + S_n*x^n/n!.
	for (n = 1; n < N; n++) {
		fterm = The(cl_LF)(fterm*x)/n;
		fsum = fsum + fterm;
		gterm = The(cl_LF)(gterm*x)/n;
		if (evenp(n))
			gterm = gterm + fterm/square((cl_I)(2*n+1));
		else
			gterm = gterm - fterm/square((cl_I)(2*n+1));
		gsum = gsum + gterm;
	}
	var cl_LF result = gsum/fsum;
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(N^2).

// Same algorithm as expintegral1, but using binary splitting to evaluate
// the sums.
const cl_LF compute_catalanconst_expintegral2 (uintC len)
{
	var uintC actuallen = len+2; // 2 guard digits
	var uintC x = (uintC)(0.693148*intDsize*actuallen)+1;
	var uintC N = (uintC)(2.718281828*x);
	CL_ALLOCA_STACK;
	var cl_pqd_series_term* args = (cl_pqd_series_term*) cl_alloca(N*sizeof(cl_pqd_series_term));
	var uintC n;
	for (n = 0; n < N; n++) {
		if (n==0) {
			init1(cl_I, args[n].p) (1);
			init1(cl_I, args[n].q) (1);
		} else {
			init1(cl_I, args[n].p) (x);
			init1(cl_I, args[n].q) (n);
		}
		init1(cl_I, args[n].d) (evenp(n)
		                        ? square((cl_I)(2*n+1))
		                        : -square((cl_I)(2*n+1)));
	}
	var cl_LF result = eval_pqd_series(N,args,actuallen);
	for (n = 0; n < N; n++) {
		args[n].p.~cl_I();
		args[n].q.~cl_I();
		args[n].d.~cl_I();
	}
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

// Using Cohen-Villegas-Zagier acceleration, but without binary splitting.
const cl_LF compute_catalanconst_cvz1 (uintC len)
{
	var uintC actuallen = len+2; // 2 guard digits
	var uintC N = (uintC)(0.39321985*intDsize*actuallen)+1;
#if 0
	var cl_LF fterm = cl_I_to_LF(2*(cl_I)N*(cl_I)N,actuallen);
	var cl_LF fsum = fterm;
	var cl_LF gterm = fterm;
	var cl_LF gsum = gterm;
	var uintC n;
	// After n loops
	//   fterm = (N+n)!N/(2n+2)!(N-n-1)!*2^(2n+2), fsum = ... + fterm,
	//   gterm = S_n*fterm, gsum = ... + gterm.
	for (n = 1; n < N; n++) {
		fterm = The(cl_LF)(fterm*(2*(cl_I)(N-n)*(cl_I)(N+n)))/((cl_I)(2*n+1)*(cl_I)(n+1));
		fsum = fsum + fterm;
		gterm = The(cl_LF)(gterm*(2*(cl_I)(N-n)*(cl_I)(N+n)))/((cl_I)(2*n+1)*(cl_I)(n+1));
		if (evenp(n))
			gterm = gterm + fterm/square((cl_I)(2*n+1));
		else
			gterm = gterm - fterm/square((cl_I)(2*n+1));
		gsum = gsum + gterm;
	}
	var cl_LF result = gsum/(cl_I_to_LF(1,actuallen)+fsum);
#else
	// Take advantage of the fact that fterm and fsum are integers.
	var cl_I fterm = 2*(cl_I)N*(cl_I)N;
	var cl_I fsum = fterm;
	var cl_LF gterm = cl_I_to_LF(fterm,actuallen);
	var cl_LF gsum = gterm;
	var uintC n;
	// After n loops
	//   fterm = (N+n)!N/(2n+2)!(N-n-1)!*2^(2n+2), fsum = ... + fterm,
	//   gterm = S_n*fterm, gsum = ... + gterm.
	for (n = 1; n < N; n++) {
		fterm = exquopos(fterm*(2*(cl_I)(N-n)*(cl_I)(N+n)),(cl_I)(2*n+1)*(cl_I)(n+1));
		fsum = fsum + fterm;
		gterm = The(cl_LF)(gterm*(2*(cl_I)(N-n)*(cl_I)(N+n)))/((cl_I)(2*n+1)*(cl_I)(n+1));
		if (evenp(n))
			gterm = gterm + cl_I_to_LF(fterm,actuallen)/square((cl_I)(2*n+1));
		else
			gterm = gterm - cl_I_to_LF(fterm,actuallen)/square((cl_I)(2*n+1));
		gsum = gsum + gterm;
	}
	var cl_LF result = gsum/cl_I_to_LF(1+fsum,actuallen);
#endif
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(N^2).

// Using Cohen-Villegas-Zagier acceleration, with binary splitting.
const cl_LF compute_catalanconst_cvz2 (uintC len)
{
	var uintC actuallen = len+2; // 2 guard digits
	var uintC N = (uintC)(0.39321985*intDsize*actuallen)+1;
	CL_ALLOCA_STACK;
	var cl_pqd_series_term* args = (cl_pqd_series_term*) cl_alloca(N*sizeof(cl_pqd_series_term));
	var uintC n;
	for (n = 0; n < N; n++) {
		init1(cl_I, args[n].p) (2*(cl_I)(N-n)*(cl_I)(N+n));
		init1(cl_I, args[n].q) ((cl_I)(2*n+1)*(cl_I)(n+1));
		init1(cl_I, args[n].d) (evenp(n)
		                        ? square((cl_I)(2*n+1))
		                        : -square((cl_I)(2*n+1)));
	}
	var cl_pqd_series_result<cl_I> sums;
	eval_pqd_series_aux(N,args,sums);
	// Here we need U/(1+S) = V/D(Q+T).
	var cl_LF result =
	  cl_I_to_LF(sums.V,actuallen) / The(cl_LF)(sums.D * cl_I_to_LF(sums.Q+sums.T,actuallen));
	for (n = 0; n < N; n++) {
		args[n].p.~cl_I();
		args[n].q.~cl_I();
		args[n].d.~cl_I();
	}
	return shorten(result,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(log(N)^2*M(N)).


const cl_LF compute_catalanconst_lupas (uintC len)
{
        // [Alexandru Lupas. Formulae for Some Classical Constants.
        //  Proceedings of ROGER-2000.
        //  http://www.lacim.uqam.ca/~plouffe/articles/alupas1.pdf]
        // [http://mathworld.wolfram.com/CatalansConstant.html]
	// G = 19/18 * sum(n=0..infty,
	//                 mul(m=1..n, -32*((80*m^3+72*m^2-18*m-19)*m^3)/
	//                             (10240*m^6+14336*m^5+2560*m^4-3072*m^3-888*m^2+72*m+27))).
	struct rational_series_stream : cl_pq_series_stream {
		cl_I n;
		static cl_pq_series_term computenext (cl_pq_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var cl_I n = thiss.n;
			var cl_pq_series_term result;
			if (zerop(n)) {
				result.p = 1;
				result.q = 1;
			} else {
				// Compute -32*((80*n^3+72*n^2-18*n-19)*n^3) (using Horner scheme):
				result.p = (19+(18+(-72-80*n)*n)*n)*n*n*n << 5;
				// Compute 10240*n^6+14336*n^5+2560*n^4-3072*n^3-888*n^2+72*n+27:
				result.q = 27+(72+(-888+(-3072+(2560+(14336+10240*n)*n)*n)*n)*n)*n;
			}
			thiss.n = plus1(n);
			return result;
		}
		rational_series_stream ()
			: cl_pq_series_stream (rational_series_stream::computenext),
			  n (0) {}
	} series;
	var uintC actuallen = len + 2; // 2 guard digits
	var uintC N = (intDsize/2)*actuallen;
	var cl_LF fsum = eval_rational_series<false>(N,series,actuallen,actuallen);
	var cl_LF g = fsum*cl_I_to_LF(19,actuallen)/cl_I_to_LF(18,actuallen);
	return shorten(g,len);
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

// Timings of the above algorithms, on an AMD Opteron 1.7 GHz, running Linux/x86.
//    N      ram    ramfast  exp1    exp2    cvz1    cvz2    lupas
//    25     0.0011  0.0010  0.0094  0.0107  0.0021  0.0016  0.0034
//    50     0.0030  0.0025  0.0280  0.0317  0.0058  0.0045  0.0095
//   100     0.0087  0.0067  0.0854  0.0941  0.0176  0.0121  0.0260
//   250     0.043   0.029   0.462   0.393   0.088   0.055   0.109
//   500     0.15    0.086   1.7     1.156   0.323   0.162   0.315
//  1000     0.57    0.25    7.5     3.23    1.27    0.487   0.864
//  2500     3.24    1.10   52.2    12.4     8.04    2.02    3.33
//  5000    13.1     3.06  218      32.7    42.1     5.59    8.80
// 10000    52.7     8.2   910      85.3   216.7    15.3    22.7
// 25000   342      29.7  6403     295    1643      54.5    77.3
// 50000            89.9                           139     195
//100000           227                             363     483
// asymp.    N^2     FAST    N^2     FAST    N^2     FAST    FAST

// (FAST means O(log(N)^2*M(N)))
//
// The "ramfast" algorithm is always fastest.

const cl_LF compute_catalanconst (uintC len)
{
	return compute_catalanconst_ramanujan_fast(len);
}
// Bit complexity (N := len): O(log(N)^2*M(N)).

const cl_LF catalanconst (uintC len)
{
	var uintC oldlen = TheLfloat(cl_LF_catalanconst())->len; // vorhandene Länge
	if (len < oldlen)
		return shorten(cl_LF_catalanconst(),len);
	if (len == oldlen)
		return cl_LF_catalanconst();

	// TheLfloat(cl_LF_catalanconst())->len um mindestens einen konstanten Faktor
	// > 1 wachsen lassen, damit es nicht zu häufig nachberechnet wird:
	var uintC newlen = len;
	oldlen += floor(oldlen,2); // oldlen * 3/2
	if (newlen < oldlen)
		newlen = oldlen;

	// gewünschte > vorhandene Länge -> muß nachberechnen:
	cl_LF_catalanconst() = compute_catalanconst(newlen);
	return (len < newlen ? shorten(cl_LF_catalanconst(),len) : cl_LF_catalanconst());
}

}  // namespace cln
