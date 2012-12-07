// square().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N square (const cl_N& x)
{
// Methode:
// x reell -> klar.
// x=a+bi -> (a^2-b^2)+(2*a*b)i
	if (realp(x)) {
		DeclareType(cl_R,x);
		return square(x);
	} else {
		DeclareType(cl_C,x);
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		return complex_C(square(a)-square(b),2*a*b);
	}
}

}  // namespace cln
