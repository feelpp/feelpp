// cl_LF_shortenwith().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "base/cl_inline2.h"
#include "float/lfloat/misc/cl_LF_precision.cc"
#include "base/cl_inline.h"
#include "float/lfloat/misc/cl_LF_exponent.cc"

namespace cln {

const cl_LF cl_LF_shortenwith (const cl_LF& x, const cl_LF& y)
{
	// Methode:
	// x = 0.0 -> Precision egal, return x.
	// ex := float_exponent(x), dx := float_digits(x), 1 ulp(x) = 2^(ex-dx).
	// ey := float_exponent(y).
	// Falls ex-dx < ey, x von Precision dx auf ex-ey verkürzen.
	var sintE ey = float_exponent_inline(y);
	var sintE ex = float_exponent_inline(x);
	var uintC dx = float_precision_inline(x);
	if (dx==0) // zerop(x) ?
		return x;
	var sintE ulpx = ex - dx;
	if ((ex<0 && ulpx>=0) // underflow?
	    || (ulpx < ey)
	   ) {	// Now ex-dx < ey, hence ex-ey < dx.
		var uintL new_dx;
		if (ex < ey)
			new_dx = intDsize*LF_minlen;
		else if ((new_dx = ex - ey) < intDsize*LF_minlen)
			new_dx = intDsize*LF_minlen;
		var uintL len = ceiling(new_dx,intDsize);
		if (intDsize*len < dx)
			return shorten(x,len);
		else
			return x;
	} else
		return x;
}

}  // namespace cln
