// zeta().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/transcendental/cl_LF_tran.h"
#include "float/lfloat/cl_LF.h"
#include "cln/integer.h"
#include "cln/exception.h"
#include "base/cl_alloca.h"

namespace cln {

const cl_LF compute_zeta_exp (int s, uintC len)
{
	// Method:
	// zeta(s) = 1/(1-2^(1-s)) sum(n=0..infty, (-1)^n/(n+1)^s),
	// with convergence acceleration through exp(x), and evaluated
	// using the binary-splitting algorithm.
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
		                        ? expt_pos(n+1,s)
		                        : -expt_pos(n+1,s));
	}
	var cl_LF result = eval_pqd_series(N,args,actuallen);
	for (n = 0; n < N; n++) {
		args[n].p.~cl_I();
		args[n].q.~cl_I();
		args[n].d.~cl_I();
	}
	result = shorten(result,len); // verkürzen und fertig
	// Zum Schluss mit 2^(s-1)/(2^(s-1)-1) multiplizieren:
	return scale_float(result,s-1) / (ash(1,s-1)-1);
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

const cl_LF compute_zeta_cvz1 (int s, uintC len)
{
	// Method:
	// zeta(s) = 1/(1-2^(1-s)) sum(n=0..infty, (-1)^n/(n+1)^s),
	// with Cohen-Villegas-Zagier convergence acceleration.
	var uintC actuallen = len+2; // 2 guard digits
	var uintC N = (uintC)(0.39321985*intDsize*actuallen)+1;
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
			gterm = gterm + cl_I_to_LF(fterm,actuallen)/expt_pos(n+1,s);
		else
			gterm = gterm - cl_I_to_LF(fterm,actuallen)/expt_pos(n+1,s);
		gsum = gsum + gterm;
	}
	var cl_LF result = gsum/cl_I_to_LF(1+fsum,actuallen);
	result = shorten(result,len); // verkürzen und fertig
	// Zum Schluss mit 2^(s-1)/(2^(s-1)-1) multiplizieren:
	return scale_float(result,s-1) / (ash(1,s-1)-1);
}
// Bit complexity (N = len): O(N^2).

const cl_LF compute_zeta_cvz2 (int s, uintC len)
{
	// Method:
	// zeta(s) = 1/(1-2^(1-s)) sum(n=0..infty, (-1)^n/(n+1)^s),
	// with Cohen-Villegas-Zagier convergence acceleration, and
	// evaluated using the binary splitting algorithm with truncation.
	var uintC actuallen = len+2; // 2 guard digits
	var uintC N = (uintC)(0.39321985*intDsize*actuallen)+1;
	struct rational_series_stream : cl_pqd_series_stream {
		uintC n;
		int s;
		uintC N;
		static cl_pqd_series_term computenext (cl_pqd_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var uintC s = thiss.s;
			var uintC N = thiss.N;
			var cl_pqd_series_term result;
			result.p = 2*(cl_I)(N-n)*(cl_I)(N+n);
			result.q = (cl_I)(2*n+1)*(cl_I)(n+1);
			result.d = evenp(n) ? expt_pos(n+1,s) : -expt_pos(n+1,s);
			thiss.n = n+1;
			return result;
		}
		rational_series_stream (int s_, uintC N_)
			: cl_pqd_series_stream (rational_series_stream::computenext),
			  n (0), s (s_), N (N_) {}
	} series(s,N);
	var cl_pqd_series_result<cl_I> sums;
	eval_pqd_series_aux(N,series,sums,actuallen);
	// Here we need U/(1+S) = V/D(Q+T).
	var cl_LF result =
	  cl_I_to_LF(sums.V,actuallen) / The(cl_LF)(sums.D * cl_I_to_LF(sums.Q+sums.T,actuallen));
	result = shorten(result,len); // verkürzen und fertig
	// Zum Schluss mit 2^(s-1)/(2^(s-1)-1) multiplizieren:
	return scale_float(result,s-1) / (ash(1,s-1)-1);
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

// Timings of the above algorithm in seconds, on a P-4, 3GHz, running Linux.
//    s                5                          15
//    N     sum_exp sum_cvz1 sum_cvz2   sum_exp sum_cvz1 sum_cvz2
//   125      0.60    0.04     0.06       1.88    0.04     0.20
//   250      1.60    0.13     0.19       4.82    0.15     0.58   
//   500      4.3     0.48     0.60      12.2     0.55     1.67
//  1000     11.0     1.87     1.63      31.7     2.11     4.60
//  2000     28.0     7.4      4.23     111       8.2     11.3
//  4000     70.2    30.6     10.6               50       44
//  8000            142       26.8              169       75
// asymp.      FAST    N^2     FAST        FAST    N^2     FAST
//
//    s               35                          75
//    N     sum_exp sum_cvz1 sum_cvz2   sum_exp sum_cvz1 sum_cvz2
//   125      4.70    0.05     0.53      11.3     0.07     1.35
//   250     12.5     0.19     1.62      28.7     0.25     3.74
//   500     31.3     0.69     4.40      70.2     0.96    10.2
//  1000     88.8     2.70    11.4      191       3.76    25.4
//  2000             10.9     28.9               15.6     64.3
//  4000             46       73                 64.4    170
//  8000            215      178                295      397
// 16000            898      419               1290      972
// asymp.      FAST    N^2     FAST        FAST    N^2     FAST
//
// The break-even point between cvz1 and cvz2 seems to grow linearly with s.

// Timings of the above algorithm, on an i486 33 MHz, running Linux.
//    s                5                           15
//    N     sum_exp sum_cvz1 sum_cvz2   sum_exp sum_cvz1 sum_cvz2
//    10       2.04    0.09    0.17        8.0     0.11    0.49
//    25       8.6     0.30    0.76       30.6     0.37    2.36
//    50      25.1     0.92    2.49       91.1     1.15    7.9
//   100               2.97    8.46                3.75   24.5
//   250              16.7    36.5                21.7   108
//   500              64.2   106                  85.3   295
//  1000             263     285                 342     788
// asymp.      FAST    N^2     FAST        FAST    N^2     FAST
//
// The break-even point between cvz1 and cvz2 seems to grow linearly with s.

const cl_LF zeta (int s, uintC len)
{
	if (!(s > 1))
		throw runtime_exception("zeta(s) with illegal s<2.");
	if (s==3)
		return zeta3(len);
	if (len < 220*(uintC)s)
		return compute_zeta_cvz1(s,len);
	else
		return compute_zeta_cvz2(s,len);
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

}  // namespace cln
