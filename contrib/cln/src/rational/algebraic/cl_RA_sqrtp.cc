// sqrtp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

bool sqrtp (const cl_RA& x, cl_RA* w)
{
// Methode:
// Bei Integers: klar.
// Bei Brüchen a/b : muß a=c^2 und b=d^2 sein. Dann ist die Wurzel = c/d
// (mit ggT(c,d)=1 und d>1).
	if (integerp(x)) {
		DeclareType(cl_I,x);
		return sqrtp(x,(cl_I*)w);
	} else {
	// x ist Ratio
	DeclareType(cl_RT,x);
	var const cl_I& b = denominator(x);
	var cl_I d;
	if (!sqrtp(b,&d)) // Nenner auf Quadratzahl testen
		return false;
	var const cl_I& a = numerator(x);
	var cl_I c;
	if (!sqrtp(a,&c)) // Zähler auf Quadratzahl testen
		return false;
	// beides Quadratzahlen -> Quotient der Wurzeln bilden
	*w = I_I_to_RT(c,d); return true;
}}

}  // namespace cln
