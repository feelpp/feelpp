// truncate1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_I truncate1 (const cl_RA& x, const cl_RA& y)
{
// Methode:
// x = a/b, y = c/d -> liefere (truncate (* a d) (* b c)).
// x Integer -> dito mit b=1.
// y Integer -> dito mit d=1.
// x und y Integer -> bekannt.
	if (integerp(x)) {
		DeclareType(cl_I,x);
		if (integerp(y)) {
			DeclareType(cl_I,y);
			return truncate1(x,y);
		} else {
			DeclareType(cl_RT,y);
			var const cl_I& c = numerator(y);
			var const cl_I& d = denominator(y);
			return truncate1(x*d,c);
		}
	} else {
		DeclareType(cl_RT,x);
		var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		if (integerp(y)) {
			DeclareType(cl_I,y);
			return truncate1(a,b*y);
		} else {
			DeclareType(cl_RT,y);
			var const cl_I& c = numerator(y);
			var const cl_I& d = denominator(y);
			return truncate1(a*d,b*c);
		}
	}
}

}  // namespace cln
