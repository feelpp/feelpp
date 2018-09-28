// cl_LF_shortenrelative().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/exception.h"

#include "base/cl_inline2.h"
#include "float/lfloat/misc/cl_LF_precision.cc"
#include "base/cl_inline.h"
#include "float/lfloat/misc/cl_LF_exponent.cc"

namespace cln {

const cl_LF cl_LF_shortenrelative (const cl_LF& x, const cl_LF& y)
{
	// Methode:
	// x = 0.0 -> Precision egal, return x.
	// ex := float_exponent(x), ey := float_exponent(y).
	// dx := float_digits(x), dy := float_digits(y).
	// 1 ulp(x) = 2^(ex-dx), 1 ulp(y) = 2^(ey-dy).
	// Falls ex-dx < ey-dy, x von Precision dx auf dy-ey+ex verkürzen.
	var sintE ey = float_exponent_inline(y);
	var sintC dy = float_precision_inline(y);
	if (dy==0) // zerop(y) ?
		throw runtime_exception();
	var sintE ex = float_exponent_inline(x);
	var sintC dx = float_precision_inline(x);
	if (dx==0) // zerop(x) ?
		return x;
	var sintE d = ex - ey;
	if (ex>=0 && ey<0 && d<0) // d overflow?
		return x;
	if (ex<0 && ey>=0 && d>=0) // d underflow?
		return LF_to_LF(x,LF_minlen);
	if (d >= dx - dy)
		return x;
	var uintC new_dx = dy + d;
	var uintC len = ceiling(new_dx,intDsize);
	if (len < LF_minlen)
		len = LF_minlen;
	if (intDsize*len < (uintC)dx)
		return shorten(x,len);
	else
		return x;
}

}  // namespace cln
