// sqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "cln/rational.h"
#include "cln/float.h"

namespace cln {

CL_INLINE const cl_R CL_INLINE_DECL(sqrt) (const cl_RA& x)
{
	var cl_RA w;
	if (sqrtp(x,&w)) // auf Quadrat testen
		return w; // war Quadrat, w ist die Wurzel
	else
		// x in Float umwandeln, dann die Wurzel ziehen:
		return sqrt(cl_float(x));
}

}  // namespace cln
