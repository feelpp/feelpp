// expt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

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
inline const cl_N expt_pos (const cl_N& x, uintL y)
{
	var cl_N a = x;
	var uintL b = y;
	while (!(b % 2)) { a = square(a); b = b >> 1; }
	var cl_N c = a;
	until (b == 1)
	  { b = b >> 1;
	    a = square(a);
	    if (b % 2) { c = a * c; }
	  }
	return c;
}

const cl_N expt (const cl_N& x, sintL y)
{
	if (realp(x)) {
		DeclareType(cl_R,x);
		// x reell -> schnellere Routine
		return expt(x,y);
	}
	if (y==0) { return 1; } // y=0 -> Ergebnis 1
	var uintL abs_y = (y<0 ? (uintL)(-y) : y); // Betrag von y nehmen
	var cl_N z = expt_pos(x,abs_y); // (expt x (abs y))
	return (y<0 ? recip(z) : z); // evtl. noch Kehrwert nehmen
}

}  // namespace cln
