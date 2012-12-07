// cos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N cos (const cl_N& x)
{
// Methode:
// x reell -> klar
// x = a+bi -> (complex (* (cos a) (cosh b)) (- (* (sin a) (sinh b))))
	if (realp(x)) {
		DeclareType(cl_R,x);
		return cos(x);
	} else {
		DeclareType(cl_C,x);
		// x=a+bi
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		var cosh_sinh_t hyp_b = cosh_sinh(b); // cosh(b), sinh(b) errechnen
		var cos_sin_t trig_a = cos_sin(a); // cos(a), sin(a) errechnen
		return complex(trig_a.cos * hyp_b.cosh, // cos(a)*cosh(b)
			       - (trig_a.sin * hyp_b.sinh) // -sin(a)*sinh(b)
			      );
	}
}

}  // namespace cln
