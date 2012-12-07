// truncate1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"

namespace cln {

const cl_I truncate1 (const cl_R& x, const cl_R& y)
{
// Methode:
// Beides rationale Zahlen -> truncate1(x,y).
// Sonst: truncate1(x/y).
	if (rationalp(x))
		if (rationalp(y)) {
			DeclareType(cl_RA,x);
			DeclareType(cl_RA,y);
			return truncate1(x,y);
		}
	return truncate1(x/y);
}

}  // namespace cln
