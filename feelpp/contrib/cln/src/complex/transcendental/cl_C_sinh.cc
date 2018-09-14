// sinh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N sinh (const cl_N& x)
{
// Methode:
// x reell -> klar
// x = a+bi -> (complex (* (sinh a) (cos b)) (* (cosh a) (sin b)))
	if (realp(x)) {
		DeclareType(cl_R,x);
		return sinh(x);
	} else {
		DeclareType(cl_C,x);
		// x=a+bi
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		var cosh_sinh_t hyp_a = cosh_sinh(a); // cosh(a), sinh(a) errechnen
		var cos_sin_t trig_b = cos_sin(b); // cos(b), sin(b) errechnen
		// Da b nicht = Fixnum 0 ist, ist auch sin(b) nicht = Fixnum 0.
		// cosh(a) /= Fixnum 0.
		return complex_C(hyp_a.sinh * trig_b.cos, // sinh(a)*cos(b)
				 hyp_a.cosh * trig_b.sin // cosh(a)*sin(b), nicht Fixnum 0
				);
	}
}

}  // namespace cln
