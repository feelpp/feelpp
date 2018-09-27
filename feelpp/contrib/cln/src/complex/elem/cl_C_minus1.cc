// minus1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N minus1 (const cl_N& x)
{
// Methode:
// x reell -> klar.
// x=a+bi -> (a-1)+bi
	if (realp(x)) {
		DeclareType(cl_R,x);
		return minus1(x);
	} else {
		DeclareType(cl_C,x);
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		return complex_C(minus1(a),b);
	}
}

}  // namespace cln
