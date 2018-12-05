// expt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_RA expt (const cl_RA& x, const cl_I& y)
{
// Methode:
// Für y>0: klar.
// Für y=0: Ergebnis 1.
// Für y<0: (/ (expt x (- y))).
	if (minusp(y))
		return recip(expt_pos(x,-y));
	elif (zerop(y))
		return 1;
	else // y > 0
		return expt_pos(x,y);
}

}  // namespace cln
