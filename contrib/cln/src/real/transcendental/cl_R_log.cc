// log().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "base/cl_N.h"
#include "real/cl_R.h"
#include "cln/rational.h"
#include "cln/float.h"

namespace cln {

const cl_R log (const cl_R& a, const cl_R& b)
{
// Methode:
// a und b rational:
//   b=1 -> Error
//   log(a,b) rational errechenbar -> liefern.
//   Sonst a und b in Floats umwandeln.
// a Float, b rational -> bei b=1 Error, sonst b := (float b a)
// a rational, b Float -> bei a=1 Ergebnis 0, sonst a := (float a b)
// a,b Floats -> log(a,b) = ln(a)/ln(b)
 {	Mutable(cl_R,a);
	Mutable(cl_R,b);
	if (rationalp(b)) {
		// b rational
		if (eq(b,1)) { throw division_by_0_exception(); }
		if (rationalp(a)) {
			// a,b beide rational
			var cl_RA l;
			if (logp(The(cl_RA)(a),The(cl_RA)(b),&l))
				return l;
			// a,b beide in Floats umwandeln:
			a = cl_float(The(cl_RA)(a)); b = cl_float(The(cl_RA)(b));
		} else
			// a Float
			b = cl_float(The(cl_RA)(b),The(cl_F)(a)); // b := (float b a)
	} else {
		// b Float
		if (rationalp(a)) {
			if (eq(a,1)) { return 0; } // a=1 -> Ergebnis 0
			a = cl_float(The(cl_RA)(a),The(cl_F)(b)); // a := (float a b)
		}
	}
	// Nun a,b beide Floats.
	return ln(The(cl_F)(a)) / ln(The(cl_F)(b));
}}

}  // namespace cln
