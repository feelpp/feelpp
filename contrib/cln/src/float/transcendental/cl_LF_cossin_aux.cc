// cl_cossin_aux().

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
#include "cln/exception.h"

#undef floor
#include <cmath>
#define floor cln_floor

namespace cln {

// Computing cos(x) = sqrt(1-sin(x)^2) instead of computing separately
// by a power series evaluation brings 20% speedup, even more for small lengths.
#define TRIVIAL_SPEEDUP

const cl_LF_cos_sin_t cl_cossin_aux (const cl_I& p, uintE lq, uintC len)
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
	// cos(p/2^lq):
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = 1, b(n) = 1,
	//   p(n) = -p^2 for n>0,
	//   q(n) = (2*n-1)*(2*n)*(2^lq)^2 for n>0.
	// sin(p/2^lq):
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = 1, b(n) = 1,
	//   p(0) = p, p(n) = -p^2 for n>0,
	//   q(0) = 2^lq, q(n) = (2*n)*(2*n+1)*(2^lq)^2 for n>0.
	var uintC actuallen = len+1; // 1 guard digit
	// How many terms do we need for M bits of precision? N/2 terms suffice,
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
	N = ceiling(N,2);
	CL_ALLOCA_STACK;
	var cl_I* pv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var cl_I* qv = (cl_I*) cl_alloca(N*sizeof(cl_I));
	var uintC n;
	var cl_I p2 = -square(p);
	var cl_LF sinsum;
	{
		init1(cl_I, pv[0]) (p);
		init1(cl_I, qv[0]) ((cl_I)1 << lq);
		for (n = 1; n < N; n++) {
			init1(cl_I, pv[n]) (p2);
			init1(cl_I, qv[n]) (((cl_I)n*(cl_I)(2*n+1)) << (2*lq+1));
		}
		var cl_pq_series series;
		series.pv = pv; series.qv = qv;
		sinsum = eval_rational_series<true>(N,series,actuallen);
		for (n = 0; n < N; n++) {
			pv[n].~cl_I();
			qv[n].~cl_I();
		}
	}
	#if !defined(TRIVIAL_SPEEDUP)
	var cl_LF cossum;
	{
		init1(cl_I, pv[0]) (1);
		init1(cl_I, qv[0]) (1);
		for (n = 1; n < N; n++) {
			init1(cl_I, pv[n]) (p2);
			init1(cl_I, qv[n]) (((cl_I)n*(cl_I)(2*n-1)) << (2*lq+1));
		}
		var cl_pq_series series;
		series.pv = pv; series.qv = qv;
		cossum = eval_rational_series<true>(N,series,actuallen);
		for (n = 0; n < N; n++) {
			pv[n].~cl_I();
			qv[n].~cl_I();
		}
	}
	#else // TRIVIAL_SPEEDUP
	var cl_LF cossum = sqrt(cl_I_to_LF(1,actuallen) - square(sinsum));
	#endif
	return cl_LF_cos_sin_t(shorten(cossum,len),shorten(sinsum,len)); // verkürzen und fertig
}}
// Bit complexity (N = len, and if p has length O(log N) and ql = O(log N)):
// O(log(N)*M(N)).

}  // namespace cln
