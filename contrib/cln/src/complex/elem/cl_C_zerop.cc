// zerop().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

bool zerop (const cl_N& x)
{
	if (realp(x)) {
		DeclareType(cl_R,x);
		return zerop(x);
	} else {
		DeclareType(cl_C,x);
		// x komplex, teste ob Real- und Imaginärteil beide = 0 sind.
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		if (zerop(a))
			if (zerop(b))
				return true;
		return false;
	}
}

}  // namespace cln
