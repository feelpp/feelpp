// round2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

const cl_RA_div_t round2 (const cl_RA& x, const cl_RA& y)
{
#if 1 // Ist das wirklich schneller??
// Methode:
// x = a/b, y = c/d ->
//   (round (* a d) (* b c)) liefert q und r.
//   Liefere q und r/(b*d).
// x Integer -> dito mit b=1.
// y Integer -> dito mit d=1.
// x und y Integer -> bekannt.
	if (integerp(x)) {
		DeclareType(cl_I,x);
		if (integerp(y)) {
			DeclareType(cl_I,y);
			var cl_I_div_t q_r = round2(x,y);
			var cl_I& q = q_r.quotient;
			var cl_I& r = q_r.remainder;
			return cl_RA_div_t(q,r);
		} else {
			DeclareType(cl_RT,y);
			var const cl_I& c = numerator(y);
			var const cl_I& d = denominator(y);
			var cl_I_div_t q_r = round2(x*d,c);
			var cl_I& q = q_r.quotient;
			var cl_I& r = q_r.remainder;
			return cl_RA_div_t(q,I_posI_div_RA(r,d));
		}
	} else {
		DeclareType(cl_RT,x);
		var const cl_I& a = numerator(x);
		var const cl_I& b = denominator(x);
		if (integerp(y)) {
			DeclareType(cl_I,y);
			var cl_I_div_t q_r = round2(a,b*y);
			var cl_I& q = q_r.quotient;
			var cl_I& r = q_r.remainder;
			return cl_RA_div_t(q,I_posI_div_RA(r,b));
		} else {
			DeclareType(cl_RT,y);
			var const cl_I& c = numerator(y);
			var const cl_I& d = denominator(y);
			var cl_I_div_t q_r = round2(a*d,b*c);
			var cl_I& q = q_r.quotient;
			var cl_I& r = q_r.remainder;
			return cl_RA_div_t(q,I_posI_div_RA(r,b*d));
		}
	}
#else
// Methode:
// round2(x/y) -> (q,r). Liefere q und x-y*q=y*r.
	var cl_RA_div_t q_r = round2(x/y);
	var cl_I& q = q_r.quotient;
	var cl_RA& r = q_r.remainder;
	return cl_RA_div_t(q,y*r);
#endif
}

}  // namespace cln
