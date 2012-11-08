// cl_N equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "base/cl_N.h"
#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

uint32 equal_hashcode (const cl_N& x)
{
	if (realp(x)) {
		DeclareType(cl_R,x);
		return equal_hashcode(x);
	} else {
		DeclareType(cl_C,x);
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		var uint32 code1 = equal_hashcode(a);
		var uint32 code2 = equal_hashcode(b);
		// Wichtig beim Kombinieren, wegen "complex canonicalization":
		// Ist imagpart=0.0, so ist der Hashcode = equal_hashcode(a).
		return code1 ^ ((code2 << 5) | (code2 >> 27));
	}
}

}  // namespace cln
