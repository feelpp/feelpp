// exp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N exp (const cl_N& x)
{
// Methode:
// x reell -> klar.
// x = a+bi -> (exp a) mit (cos b) + i (sin b) multiplizieren:
//             (complex (* (exp a) (cos b)) (* (exp a) (sin b)))
	if (realp(x)) {
		DeclareType(cl_R,x);
		return exp(x);
	} else {
		DeclareType(cl_C,x);
		// x=a+bi
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		var cos_sin_t unitvec = cos_sin(b); // (cos b), (sin b)
		// Da b nicht = Fixnum 0 ist, ist auch sin(b) nicht = Fixnum 0.
		var cl_R exp_a = exp(a); // (exp a)
		return complex_C(exp_a * unitvec.cos, exp_a * unitvec.sin);
	}
}

}  // namespace cln
