// expt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

// Methode:
// Für y>0:
//   a:=x, b:=y.
//   Solange b gerade, setze a:=a*a, b:=b/2. [a^b bleibt invariant, = x^y.]
//   c:=a.
//   Solange b:=floor(b/2) >0 ist,
//     setze a:=a*a, und falls b ungerade, setze c:=a*c.
//   Ergebnis c.
// Für y=0: Ergebnis 1.
// Für y<0: (/ (expt x (- y))).

// Assume y>0.
inline const cl_R expt_pos (const cl_R& x, const cl_I& y)
{
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		return expt(x,y); // x rational -> schnellere Routine
	} else {
		DeclareType(cl_F,x);
		var cl_F a = x;
		var cl_I b = y;
		while (!oddp(b)) { a = square(a); b = b >> 1; }
		var cl_F c = a;
		until (eq(b,1))
		  { b = b >> 1;
		    a = square(a);
		    if (oddp(b)) { c = a * c; }
		  }
		return c;
	}
}

const cl_R expt (const cl_R& x, const cl_I& y)
{
	if (eq(y,0)) { return 1; } // y=0 -> Ergebnis 1
	var bool y_negative = minusp(y);
	var cl_I abs_y = (y_negative ? -y : y); // Betrag von y nehmen
	var cl_R z = expt_pos(x,abs_y); // (expt x (abs y))
	return (y_negative ? recip(z) : z); // evtl. noch Kehrwert nehmen
}

}  // namespace cln
