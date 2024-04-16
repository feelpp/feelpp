// zeta3().

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

const cl_LF zeta3 (uintC len)
{
	struct rational_series_stream : cl_pqa_series_stream {
		uintC n;
		static cl_pqa_series_term computenext (cl_pqa_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_pqa_series_term result;
			if (n==0) {
				result.p = 1;
			} else {
				result.p = -expt_pos(n,5);
			}
			result.q = expt_pos(2*n+1,5)<<5;
			result.a = 205*square((cl_I)n) + 250*(cl_I)n + 77;
			thiss.n = n+1;
			return result;
		}
		rational_series_stream ()
			: cl_pqa_series_stream (rational_series_stream::computenext),
			  n (0) {}
	} series;
	// Method:
	//            /infinity                                  \ 
	//            | -----       (n + 1)       2              |
	//        1   |  \      (-1)        (205 n  - 160 n + 32)|
	//        -   |   )     ---------------------------------|
	//        2   |  /              5                 5      |
	//            | -----          n  binomial(2 n, n)       |
	//            \ n = 1                                    /
	//
	// The formula used to compute Zeta(3) has reference in the paper
	// "Hypergeometric Series Acceleration via the WZ method" by
	// T. Amdeberhan and Doron Zeilberger,
	// Electronic J. Combin. 4 (1997), R3.
	//
	// Computation of the sum:
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = 205*n^2+250*n+77, b(n) = 1,
	//   p(0) = 1, p(n) = -n^5 for n>0, q(n) = 32*(2n+1)^5.
	var uintC actuallen = len+2; // 2 guard digits
	var uintC N = ceiling(actuallen*intDsize,10);
	// 1024^-N <= 2^(-intDsize*actuallen).
	var cl_LF sum = eval_rational_series<false>(N,series,actuallen,actuallen);
	return scale_float(shorten(sum,len),-1);
}
// Bit complexity (N := len): O(log(N)^2*M(N)).

// Timings of the above algorithm, on an i486 33 MHz, running Linux.
//    N   sum_exp sum_cvz1 sum_cvz2 hypgeom
//    10     1.17    0.081   0.125   0.013
//    25     5.1     0.23    0.50    0.045
//    50    15.7     0.66    1.62    0.14
//   100    45.5     1.93    5.4     0.44
//   250   169      13.1    25.1     2.03
//   500   436      56.5    70.6     6.44
//  1000           236     192      18.2
//  2500                            78.3
//  5000                           202
// 10000                           522
// 25000                          1512
// 50000                          3723
// asymp.    FAST     N^2    FAST    FAST
// (FAST means O(log(N)^2*M(N)))

}  // namespace cln
