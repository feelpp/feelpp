// phase().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_R phase (const cl_N& x)
{
// Methode:
// (= x 0) -> willkürliches Ergebnis 0
// x reell -> Winkel von (x,0) in Polarkoordinaten
// x komplex -> Winkel von ((realpart x),(imagpart x)) in Polarkoordinaten
	if (zerop(x))
		return 0;
	if (realp(x)) {
		DeclareType(cl_R,x);
		return atan(x,0);
	} else {
		DeclareType(cl_C,x);
		return atan(realpart(x),imagpart(x));
	}
}

}  // namespace cln
