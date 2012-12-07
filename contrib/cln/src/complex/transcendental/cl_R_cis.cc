// cis().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N cis (const cl_R& x)
{
// Methode:
// (complex (cos x) (sin x))
	var cos_sin_t trig = cos_sin(x);
	return complex(trig.cos, trig.sin);
}

}  // namespace cln
