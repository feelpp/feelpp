// cl_atan_recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/integer.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"
#include "float/transcendental/cl_LF_tran.h"

#undef floor
#include <cmath>
#define floor cln_floor

namespace cln {

// Method:
// See examples/atan_recip.cc for a comparison of the algorithms.
// Here we take algorithm 2d. It's the fastest throughout the range.

const cl_LF cl_atan_recip (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var cl_I m2 = m*m+1;
	var uintC N = (uintC)(0.69314718*intDsize*actuallen/::log(double_approx(m2))) + 1;
	struct rational_series_stream : cl_pq_series_stream {
		var uintC n;
		var cl_I m;
		var cl_I m2;
		static cl_pq_series_term computenext (cl_pq_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_pq_series_term result;
			if (n==0) {
				result.p = thiss.m;
				result.q = thiss.m2;
			} else {
				result.p = 2*n;
				result.q = (2*n+1)*thiss.m2;
			}
			thiss.n = n+1;
			return result;
		}
		rational_series_stream(const cl_I& m_, const cl_I& m2_)
			: cl_pq_series_stream (rational_series_stream::computenext),
			  n(0), m(m_), m2(m2_) {}
	} series(m,m2);
	var cl_LF result = eval_rational_series<false>(N,series,actuallen);
	return shorten(result,len);
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

}  // namespace cln
