// exp1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/transcendental/cl_LF_tran.h"
#include "float/lfloat/cl_LF.h"
#include "cln/integer.h"

#undef floor
#include <cmath>
#define floor cln_floor

namespace cln {

const cl_LF compute_exp1 (uintC len)
{
	// Evaluate a sum(0 <= n < N, a(n)/b(n) * (p(0)...p(n))/(q(0)...q(n)))
	// with appropriate N, and
	//   a(n) = 1, b(n) = 1, p(n) = 1, q(n) = n for n>0.
	var uintC actuallen = len+1; // 1 guard digit
	// How many terms do we need for M bits of precision? N terms suffice,
	// provided that
	//   1/N! < 2^-M
	// <==   N*(log(N)-1) > M*log(2)
	// First approximation:
	//   N0 = M will suffice, so put N<=N0.
	// Second approximation:
	//   N1 = floor(M*log(2)/(log(N0)-1)), slightly too small, so put N>=N1.
	// Third approximation:
	//   N2 = ceiling(M*log(2)/(log(N1)-1)), slightly too large.
	//   N = N2+2, two more terms for safety.
	var uintC N0 = intDsize*actuallen;
	var uintC N1 = (uintC)(0.693147*intDsize*actuallen/(::log((double)N0)-1.0));
	var uintC N2 = (uintC)(0.693148*intDsize*actuallen/(::log((double)N1)-1.0))+1;
	var uintC N = N2+2;
	struct rational_series_stream : cl_q_series_stream {
		var uintC n;
		static cl_q_series_term computenext (cl_q_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_q_series_term result;
			result.q = (n==0 ? 1 : n);
			thiss.n = n+1;
			return result;
		}
		rational_series_stream()
			: cl_q_series_stream (rational_series_stream::computenext),
			  n(0) {}
	} series;
	var cl_LF fsum = eval_rational_series<false>(N,series,actuallen);
	return shorten(fsum,len); // verkürzen und fertig
}
// Bit complexity (N = len): O(log(N)*M(N)).

// Timings of the above algorithm, on an i486 33 MHz, running Linux.
//    N      
//    10     0.0040
//    25     0.0096
//    50     0.0218
//   100     0.057
//   250     0.24
//   500     0.78
//  1000     2.22
//  2500     7.6
//  5000    17.8
// 10000    41.4
// 25000   111
// 50000   254

const cl_LF exp1 (uintC len)
{
	var uintC oldlen = TheLfloat(cl_LF_exp1())->len; // vorhandene Länge
	if (len < oldlen)
		return shorten(cl_LF_exp1(),len);
	if (len == oldlen)
		return cl_LF_exp1();

	// TheLfloat(cl_LF_exp1())->len um mindestens einen konstanten Faktor
	// > 1 wachsen lassen, damit es nicht zu häufig nachberechnet wird:
	var uintC newlen = len;
	oldlen += floor(oldlen,2); // oldlen * 3/2
	if (newlen < oldlen)
		newlen = oldlen;

	// gewünschte > vorhandene Länge -> muß nachberechnen:
	cl_LF_exp1() = compute_exp1(newlen); // (exp 1)
	return (len < newlen ? shorten(cl_LF_exp1(),len) : cl_LF_exp1());
}

}  // namespace cln
