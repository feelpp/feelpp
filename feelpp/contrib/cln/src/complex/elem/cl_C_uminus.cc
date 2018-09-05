// unary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N operator- (const cl_N& x)
{
// Methode:
// x reell -> klar.
// x=a+bi -> (-a) + (-b) i
	if (realp(x)) {
		DeclareType(cl_R,x);
		return -x;
	} else {
		DeclareType(cl_C,x);
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		return complex_C(-a,-b);
	}
}

}  // namespace cln
