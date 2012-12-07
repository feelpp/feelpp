// round1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"

namespace cln {

const cl_I round1 (const cl_R& x, const cl_R& y)
{
// Methode:
// Beides rationale Zahlen -> round1(x,y).
// Sonst: round1(x/y).
	if (rationalp(x))
		if (rationalp(y)) {
			DeclareType(cl_RA,x);
			DeclareType(cl_RA,y);
			return round1(x,y);
		}
	return round1(x/y);
}

}  // namespace cln
