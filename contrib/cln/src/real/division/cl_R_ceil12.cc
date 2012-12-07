// ceiling1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"

namespace cln {

const cl_I ceiling1 (const cl_R& x, const cl_R& y)
{
// Methode:
// Beides rationale Zahlen -> ceiling1(x,y).
// Sonst: ceiling1(x/y).
	if (rationalp(x))
		if (rationalp(y)) {
			DeclareType(cl_RA,x);
			DeclareType(cl_RA,y);
			return ceiling1(x,y);
		}
	return ceiling1(x/y);
}

}  // namespace cln
