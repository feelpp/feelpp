// I_posI_div_RA().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "rational/cl_RA.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_RA I_posI_div_RA (const cl_I& a, const cl_I& b)
{
// Methode:
// d:=ggT(a,b).
// Falls d=1: I_I_to_RA anwenden,
// sonst: I_I_to_RA auf a/d und b/d anwenden.
	var cl_I d = gcd(a,b);
	if (eq(d,1))
		return I_I_to_RA(a,b);
	else
		return I_I_to_RA(exquo(a,d),exquopos(b,d));
}

}  // namespace cln
