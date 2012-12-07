// cl_RA_LF_div().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/lfloat.h"
#include "rational/cl_RA.h"

namespace cln {

const cl_R cl_RA_LF_div (const cl_RA& x, const cl_LF& y)
{
// Method:
// Write x = u/v. Return u/(v*y).
	if (integerp(x)) {
		DeclareType(cl_I,x);
		return cl_I_LF_div(x,y);
	} else {
		DeclareType(cl_RT,x);
		var const cl_I& u = TheRatio(x)->numerator;
		var const cl_I& v = TheRatio(x)->denominator; // v /= 0
		return cl_I_LF_div(u,The(cl_LF)(cl_LF_I_mul(y,v)));
	}
}

// Timings on an i486 33 MHz, running Linux.
// First timing:  u/(v*y), using cl_LF_I_mul and cl_LF_I_div
// Second timing: (u/v)/y, using cl_RA_to_LF and operator/
// With x_length = 100, in 0.01 sec.
//     num_length    50          70          100         200           500
// den_length
//
//         50     0.82 0.98   1.38 1.61   2.48 2.80   8.37  9.05   46.92 48.59
//
//         70     0.81 1.16   1.43 1.85   2.65 3.18   8.64  9.76   47.53 50.34
//
//        100     0.82 1.43   1.44 2.24   2.79 3.69   8.98 10.77   48.42 52.81
//
//        200     0.80 2.31   1.43 3.44   2.78 5.43   9.93 14.31   51.08 61.57
//
//        500     0.82 4.62   1.44 6.21   2.76 9.40   9.97 22.66   55.52 87.71
//
// With x_length = 1000, in sec.
//     num_length    500         700         1000        2000         5000
// den_length
//
//        500     0.55 0.88   0.96 1.39   1.60 2.24   4.35 5.60   10.83 14.18
//
//        700     0.55 0.95   0.98 1.57   1.65 2.58   4.47 6.34   10.61 15.88
//
//       1000     0.56 1.00   0.98 1.67   1.69 2.70   4.61 7.13   10.65 16.29
//
//       2000     0.56 1.25   0.98 1.97   1.69 3.05   4.88 7.74   10.75 17.03
//
//       5000     0.55 1.94   0.98 2.30   1.70 3.33   4.90 7.74   11.94 19.29
//
// We see that the first approach is always better than the second.

}  // namespace cln
