// unary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA operator- (const cl_RA& r)
{
// Methode:
// r Integer -> klar.
// r = a/b -> Ergebnis (- a)/b
	if (integerp(r)) {
		DeclareType(cl_I,r);
		return -r;
	} else {
		DeclareType(cl_RT,r);
		var const cl_I& a = numerator(r);
		var const cl_I& b = denominator(r);
		// Immer noch b>1 und ggT(-a,b) = ggT(a,b) = 1
		return I_I_to_RT(-a,b);
	}
}

}  // namespace cln
