// equal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

bool equal (const cl_RA& r, const cl_RA& s)
{
// Methode:
// r,s Integer -> klar
// r Integer, s Ratio: verschieden.
// r Ratio, s Integer: verschieden.
// r,s Ratios: r=a/b, s=c/d. Vergleiche a mit c und b mit d.
	if (integerp(r))
		if (integerp(s)) {
			// r,s beide Integer
			DeclareType(cl_I,r);
			DeclareType(cl_I,s);
			return equal(r,s);
		} else
			// r Integer, s Ratio
			return false;
	else
		if (integerp(s))
			// r Ratio, s Integer
			return false;
		else {
			DeclareType(cl_RT,r);
			DeclareType(cl_RT,s);
			// r,s Ratios
			if (!equal(numerator(r),numerator(s)))
				return false;
			if (!equal(denominator(r),denominator(s)))
				return false;
			return true;
		}
}

}  // namespace cln
