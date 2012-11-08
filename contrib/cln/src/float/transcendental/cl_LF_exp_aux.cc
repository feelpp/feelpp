// cl_exp_aux().

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

#undef floor
#include <cmath>
#define floor cln_floor

namespace cln {

const cl_LF cl_exp_aux (const cl_I& p, uintE lq, uintC len)
{
 {	Mutable(cl_I,p);
	var uintE lp = integer_length(p); // now |p| < 2^lp.
	if (!(lp <= lq)) throw runtime_exception();
	lp = lq - lp; // now |p/2^lq| < 2^-lp.
	// Minimize lq (saves computation time).
	{
		var uintC lp2 = ord2(p);
		if (lp2 > 0) {
			p = p >> lp2;
			lq = lq - lp2;
		}
	}
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = 1, b(n) = 1, p(n) = p for n>0, q(n) = n*2^lq for n>0.
	var uintC actuallen = len+1; // 1 guard digit
	// How many terms do we need for M bits of precision? N terms suffice,
	// provided that
	//   1/(2^(N*lp)*N!) < 2^-M
	// <==   N*(log(N)-1)+N*lp*log(2) > M*log(2)
	// First approximation:
	//   N0 = M will suffice, so put N<=N0.
	// Second approximation:
	//   N1 = floor(M*log(2)/(log(N0)-1+lp*log(2))), slightly too small,
	//   so put N>=N1.
	// Third approximation:
	//   N2 = ceiling(M*log(2)/(log(N1)-1+lp*log(2))), slightly too large.
	//   N = N2+2, two more terms for safety.
	var uintC N0 = intDsize*actuallen;
	var uintC N1 = (uintC)(0.693147*intDsize*actuallen/(::log((double)N0)-1.0+0.693148*lp));
	var uintC N2 = (uintC)(0.693148*intDsize*actuallen/(::log((double)N1)-1.0+0.693147*lp))+1;
	var uintC N = N2+2;
	struct rational_series_stream : cl_pq_series_stream {
		var uintC n;
		var cl_I p;
		var uintE lq;
		static cl_pq_series_term computenext (cl_pq_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_pq_series_term result;
			if (n==0) {
				result.p = 1;
				result.q = 1;
			} else {
				result.p = thiss.p;
				result.q = (cl_I)n << thiss.lq;
			}
			thiss.n = n+1;
			return result;
		}
		rational_series_stream(const cl_I& p_, uintE lq_)
			: cl_pq_series_stream (rational_series_stream::computenext),
			  n (0), p(p_), lq(lq_) {}
	} series(p, lq);
	var cl_LF fsum = eval_rational_series<true>(N,series,actuallen);
	return shorten(fsum,len); // verkürzen und fertig
}}
// Bit complexity (N = len, and if p has length O(log N) and ql = O(log N)):
// O(log(N)*M(N)).

}  // namespace cln
