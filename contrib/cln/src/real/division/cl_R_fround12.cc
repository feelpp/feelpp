// fround().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"

namespace cln {

const cl_F fround (const cl_R& x, const cl_R& y)
{
// Methode:
// x,y beide rational: round(x,y), Quotienten in Float umwandeln.
// Sonst: fround(x/y).
	if (rationalp(x))
		if (rationalp(y)) {
			// beides rationale Zahlen
			DeclareType(cl_RA,x);
			DeclareType(cl_RA,y);
			return cl_float(round1(x,y));
		}
	return fround(x/y);
}

}  // namespace cln
