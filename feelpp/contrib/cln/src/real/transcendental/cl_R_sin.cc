// sin().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "cln/float.h"
#include "real/cl_R.h"

namespace cln {

const cl_R sin (const cl_R& x)
{
// Methode:
// x rational -> bei x=0 0 als Ergebnis, sonst x in Float umwandeln.
// x Float -> bekannt.

	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		if (zerop(x)) // x=0 -> 0 als Ergebnis
			return 0;
		return sin(cl_float(x)); // sonst in Float umwandeln
	} else {
		DeclareType(cl_F,x);
		return sin(x);
	}
}

}  // namespace cln
