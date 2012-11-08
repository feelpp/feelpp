// cis().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N cis (const cl_N& x)
{
// Methode:
// x reell -> (complex (cos x) (sin x))
// x = a+bi -> (complex (* (exp (- b)) (cos a)) (* (exp (- b)) (sin a)))
	if (realp(x)) {
		DeclareType(cl_R,x);
		var cos_sin_t trig = cos_sin(x);
		return complex(trig.cos, trig.sin);
	} else {
		DeclareType(cl_C,x);
		// x=a+bi
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		var cos_sin_t trig_a = cos_sin(a); // cos(a), sin(a) errechnen
		var cl_R exp_minusb = exp(-b); // (exp (- b))
		return complex(exp_minusb*trig_a.cos, // (* (exp (- b)) (cos a))
			       exp_minusb*trig_a.sin); // (* (exp (- b)) (sin a))
	}
}

}  // namespace cln
