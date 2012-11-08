// cl_atanh_recip().

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
// See examples/atanh_recip.cc for a comparison of the algorithms.
// Here we take algorithm 1d. It's the fastest throughout the range.

const cl_LF cl_atanh_recip (cl_I m, uintC len)
{
	var uintC actuallen = len + 1;
	var uintC N = (uintC)(0.69314718*intDsize/2*actuallen/::log(double_approx(m))) + 1;
	struct rational_series_stream : cl_qb_series_stream {
		var uintC n;
		var cl_I m;
		var cl_I m2;
		static cl_qb_series_term computenext (cl_qb_series_stream& thisss)
		{
			var rational_series_stream& thiss = (rational_series_stream&)thisss;
			var uintC n = thiss.n;
			var cl_qb_series_term result;
			result.b = 2*n+1;
			result.q = (n==0 ? thiss.m : thiss.m2);
			thiss.n = n+1;
			return result;
		}
		rational_series_stream(const cl_I& m_)
			: cl_qb_series_stream (rational_series_stream::computenext),
			  n(0), m(m_), m2(square(m_)) {}
	} series(m);
	var cl_LF result = eval_rational_series<false>(N,series,actuallen);
	return shorten(result,len);
}
// Bit complexity (N = len): O(log(N)^2*M(N)).

}  // namespace cln
