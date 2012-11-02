// square().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA square (const cl_RA& r)
{
// Methode:
// r Integer -> klar.
// r = a/b -> Ergebnis a^2/b^2
	if (integerp(r)) {
		DeclareType(cl_I,r);
		return square(r);
	} else {
		DeclareType(cl_RT,r);
		var const cl_I& a = numerator(r);
		var const cl_I& b = denominator(r);
		// Immer noch b^2>1 und ggT(a^2,b^2) = ggT(a,b)^2 = 1
		return I_I_to_RT(square(a),square(b));
	}
}

}  // namespace cln
