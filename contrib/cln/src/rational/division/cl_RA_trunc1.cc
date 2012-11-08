// truncate1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_I truncate1 (const cl_RA& x)
{
// Methode:
// x Integer -> q := x
// x Ratio a/b -> (truncate a b)
	if (integerp(x)) {
		DeclareType(cl_I,x);
		return x;
	} else {
		DeclareType(cl_RT,x);
		var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		return truncate1(a,b);
	}
}

}  // namespace cln
