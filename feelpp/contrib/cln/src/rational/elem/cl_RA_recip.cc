// recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "cln/exception.h"
#include "rational/cl_RA.h"
#include "integer/cl_I.h"

namespace cln {

const cl_RA recip (const cl_RA& r)
{
// Methode:
// r=0 -> Error.
// a:=(numerator r), b:=(denominator r).
// a>0 -> Ergebnis b/a (mit ggT(b,a)=1).
// a<0 -> Ergebnis (- b)/(- a) (mit ggT(-b,-a)=1).
	if (zerop(r))
		throw division_by_0_exception();
	var cl_I a;
	var cl_I b;
	RA_numden_I_I(r, a = , b = );
	if (!minusp(a))
		return I_I_to_RA(b,a);
	else
		// a<0
		return I_I_to_RA(-b,-a);
}

}  // namespace cln
