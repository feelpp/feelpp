// cosh_sinh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "cln/float.h"
#include "real/cl_R.h"

namespace cln {

const cosh_sinh_t cosh_sinh (const cl_R& x)
{
// Methode:
// x rational -> bei x=0 (1,0) als Ergebnis, sonst x in Float umwandeln.
// x Float -> bekannt.

	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		if (zerop(x)) // x=0 -> (1,0) als Ergebnis
			return cosh_sinh_t(1,0);
		return cosh_sinh(cl_float(x)); // sonst in Float umwandeln
	} else {
		DeclareType(cl_F,x);
		return cosh_sinh(x);
	}
}

}  // namespace cln
