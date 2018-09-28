// cl_RA equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "base/cl_N.h"
#include "rational/cl_RA.h"
#include "base/cl_inline.h"
#include "integer/misc/cl_I_eqhashcode.cc"

namespace cln {

CL_INLINE2 uint32 CL_INLINE2_DECL(equal_hashcode) (const cl_RA& r)
{
	if (integerp(r)) {
		DeclareType(cl_I,r);
		return equal_hashcode_inline(r);
	} else {
		// Making sure that a float and its rational equivalent have
		// the same hash code is tricky. This code depends on the fact
		// that the equal_hashcode_low macro is linear in `exp'.
		DeclareType(cl_RT,r);
		return equal_hashcode_inline(numerator(r))
		       - equal_hashcode_inline(denominator(r))
		       + equal_hashcode_one;
	}
}

}  // namespace cln
