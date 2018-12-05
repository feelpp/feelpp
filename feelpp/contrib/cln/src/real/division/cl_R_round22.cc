// round2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "real/division/cl_R_div_t.h"

namespace cln {

const cl_R_div_t round2 (const cl_R& x, const cl_R& y)
{
// Methode:
// Beides rationale Zahlen -> round2(x,y).
// Sonst: round2(x/y) -> (q,r). Liefere q und x-y*q=y*r.
	if (rationalp(x))
		if (rationalp(y)) {
			DeclareType(cl_RA,x);
			DeclareType(cl_RA,y);
			return round2(x,y);
		}
	var cl_R_div_t q_r = round2(x/y);
	var cl_I& q = q_r.quotient;
	var cl_R& r = q_r.remainder;
	return cl_R_div_t(q,y*r);
}

}  // namespace cln
