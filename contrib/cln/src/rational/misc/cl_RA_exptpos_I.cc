// expt_pos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA expt_pos (const cl_RA& x, const cl_I& y)
{
  // Methode:
  // x Integer -> klar
  // x Ratio a/b -> x^y = (a^y)/(b^y), gekürzt, mit b^y>=b>1.
	if (integerp(x)) {
		DeclareType(cl_I,x);
		return expt_pos(x,y);
	} else {
		DeclareType(cl_RT,x);
		var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		return I_I_to_RT(expt_pos(a,y),expt_pos(b,y));
	}
}

}  // namespace cln
