// fceiling2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "real/division/cl_R_div_t.h"

namespace cln {

const cl_R_fdiv_t fceiling2 (const cl_R& x, const cl_R& y)
{
// Methode:
// x,y beide rational: ceiling(x,y), Quotienten in Float umwandeln.
// Sonst: fceiling(x/y) -> q,r. Liefere die Werte q und x-y*q = y*r.
	if (rationalp(x))
		if (rationalp(y)) {
			// beides rationale Zahlen
			DeclareType(cl_RA,x);
			DeclareType(cl_RA,y);
			var cl_R_div_t q_r = ceiling2(x,y);
			var cl_I& q = q_r.quotient;
			var cl_R& r = q_r.remainder;
			return cl_R_fdiv_t(cl_float(q),r);
		}
	var cl_R_fdiv_t q_r = fceiling2(x/y);
	var cl_F& q = q_r.quotient;
	var cl_R& r = q_r.remainder;
	return cl_R_fdiv_t(q,y*r);
}

}  // namespace cln
