// binary operator /

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA operator/ (const cl_RA& r, const cl_RA& s)
{
// Methode:
// (* r (/ s))
	if (integerp(r) && integerp(s)) {
		DeclareType(cl_I,r);
		DeclareType(cl_I,s);
		// r und s Integers -> schnell abhandeln
		return I_I_div_RA(r,s);
	}
	return r * recip(s);
}

}  // namespace cln
