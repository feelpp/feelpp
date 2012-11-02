// rem().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/integer.h"

namespace cln {

const cl_R rem (const cl_R& x, const cl_R& y)
{
// Methode:
// Beides Integers -> rem(x,y).
// Sonst: truncate2(x/y) -> (q,r). Liefere x-y*q=y*r.
	if (integerp(x))
		if (integerp(y)) {
			DeclareType(cl_I,x);
			DeclareType(cl_I,y);
			return rem(x,y);
		}
	return y * truncate2(x/y).remainder;
}

}  // namespace cln
