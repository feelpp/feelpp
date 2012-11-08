// cl_LF_RA_mul().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/lfloat.h"
#include "rational/cl_RA.h"

namespace cln {

const cl_R cl_LF_RA_mul (const cl_LF& x, const cl_RA& y)
{
// Method:
// Write y = u/v. Return (x*u)/v.
	if (integerp(y)) {
		DeclareType(cl_I,y);
		return cl_LF_I_mul(x,y);
	} else {
		DeclareType(cl_RT,y);
		var const cl_I& u = TheRatio(y)->numerator; // u /= 0
		var const cl_I& v = TheRatio(y)->denominator; // v /= 0
		return cl_LF_I_div(The(cl_LF)(cl_LF_I_mul(x,u)),v);
	}
}

// Timings on an i486 33 MHz, running Linux, in 0.01 sec.
// First timing:  (x*u)/v, using cl_LF_I_mul and cl_LF_I_div
// Second timing: x*(u/v), using cl_RA_to_LF and operator*
// with x_length = 100.
//     num_length    50          70          100         200         500
// den_length
//
//         50     1.59 1.92   1.76 1.92   1.89 1.92   1.90 2.78   1.93 8.07
//
//         70     1.93 2.26   2.14 2.25   2.25 2.25   2.27 2.78   2.25 8.07
//
//        100     2.44 2.77   2.65 2.77   2.79 2.80   2.80 2.81   2.80 8.05
//
//        200     2.46 4.53   2.65 4.54   2.78 4.55   2.79 4.55   2.77 8.05
//
//        500     2.45 8.38   2.65 8.50   2.76 8.49   2.79 8.52   2.81 8.55
//
// We see that the first approach is always better than the second, except if
//    den_length = x_length && x_length <= num_length <= 2*x_length
// when both are equally fast.

}  // namespace cln
