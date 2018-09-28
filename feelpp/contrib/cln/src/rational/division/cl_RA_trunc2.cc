// truncate2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA_div_t truncate2 (const cl_RA& x)
{
// Methode:
// x Integer -> (q,r) := (x,0)
// x Ratio a/b ->
//   (truncate a b) liefert q und r.
//   Liefere q und r/b (mit b>1 und ggT(r,b)=ggT(r+q*b,b)=ggT(a,b)=1).
	if (integerp(x)) {
		DeclareType(cl_I,x);
		// (q,r) := (x,0)
		return cl_RA_div_t(x,0);
	} else {
		DeclareType(cl_RT,x);
		var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		var cl_I_div_t q_r = truncate2(a,b);
		var cl_I& q = q_r.quotient;
		var cl_I& r = q_r.remainder;
		return cl_RA_div_t(q,I_I_to_RT(r,b));
	}
}

}  // namespace cln
