// cl_F_shortenrelative().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "cln/exception.h"

namespace cln {

const cl_F cl_F_shortenrelative (const cl_F& x, const cl_F& y)
{
	// Methode:
	// x = 0.0 -> Precision egal, return x.
	// ex := float_exponent(x), ey := float_exponent(y).
	// dx := float_digits(x), dy := float_digits(y).
	// 1 ulp(x) = 2^(ex-dx), 1 ulp(y) = 2^(ey-dy).
	// Falls ex-dx < ey-dy, x von Precision dx auf dy-ey+ex verkürzen.
	var sintE ey = float_exponent(y);
	var sintC dy = float_precision(y);
	if (dy==0) // zerop(y) ?
		throw runtime_exception();
	var sintE ex = float_exponent(x);
	var sintC dx = float_precision(x);
	if (dx==0) // zerop(x) ?
		return x;
	var sintE d = ex - ey;
	if (ex>=0 && ey<0 && d<0) // d overflow?
		return x;
	if (ex<0 && ey>=0 && d>=0) // d underflow?
		return cl_F_to_SF(x);
	if (d >= dx - dy)
		return x;
	var uintC new_dx = dy + d;
	floatformatcase(new_dx
	,	return cl_F_to_SF(x);
	,	return cl_F_to_FF(x);
	,	return cl_F_to_DF(x);
	,	if (intDsize*len < (uintC)dx)
			return shorten(The(cl_LF)(x),len);
		else
			return x;
	);
}

}  // namespace cln
