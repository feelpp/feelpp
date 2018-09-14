// logp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "integer/cl_I.h"
#include "rational/cl_RA.h"

namespace cln {

bool logp (const cl_RA& a, const cl_RA& b, cl_RA* pl)
{
// Methode:
// a=1 -> Ergebnis 0
// b Integer:
//   a Integer: log(a,b) rational errechenbar -> liefern
//   a Ratio: a=a1/a2 mit a1>0, a2>1.
//            a1=1 und log(a2,b) rational errechenbar -> -log(a2,b) liefern
// b Ratio: a=a1/a2, b=b1/b2 mit a1>0, a2>0, b1>0, b2>1.
//          log(a2,b2) rational errechenbar ->
//             b1=1 -> bei a1=1 liefern, sonst nicht.
//             b1>1 -> log(a1,b1) rational errechenbar und
//                     log(a1,b1)=log(a2,b2) -> liefern, sonst nicht.
//          sonst a1,a2 vertauschen:
//            log(a2/a1,b1/b2) versuchen (wie oben) ->
//              -log(a2/a1,b1/b2) liefern

	if (eq(a,1)) { // a=1 -> Ergebnis 0
		*pl = 0; return true;
	}
	if (integerp(b)) {
		// b Integer
		DeclareType(cl_I,b);
		if (integerp(a)) {
			// a,b beide Integers
			DeclareType(cl_I,a);
			return logp(a,b,pl);
		} else {
			// a Ratio, b Integer
			DeclareType(cl_RT,a);
			var const cl_I& a1 = numerator(a);
			var const cl_I& a2 = denominator(a);
			if (!eq(a1,1))
				return false;
			// a1=1
			var cl_RA l;
			if (logp(a2,b,&l)) {
				*pl = -l; return true;
			} else
				return false;
		}
	} else {
		// a rational, b Ratio
		DeclareType(cl_RT,b);
		var cl_I a1;
		var cl_I a2;
		RA_numden_I_I(a, a1 =, a2 =);
		var const cl_I& b1 = numerator(b);
		var const cl_I& b2 = denominator(b);
		{
			var cl_RA l2;
			// rationalen log(a2,b2) versuchen
			if (logp(a2,b2,&l2)) {
				if (eq(b1,1)) {
					if (eq(a1,1))
						{ *pl = l2; return true; }
					else
						return false;
				} else {
					var cl_RA l1;
					// rationalen log(a1,b1) versuchen
					if (logp(a1,b1,&l1))
						if (l1 == l2)
							{ *pl = l2; return true; }
					return false;
				}
			}
		}
		{
			var cl_RA l2;
			// rationalen log(a1,b2) versuchen
			if (logp(a1,b2,&l2)) {
				if (eq(b1,1)) {
					if (eq(a2,1))
						{ *pl = -l2; return true; }
					else
						return false;
				} else {
					var cl_RA l1;
					// rationalen log(a2,b1) versuchen
					if (logp(a2,b1,&l1))
						if (l1 == l2)
							{ *pl = -l2; return true; }
				}
			}
		}
		return false;
	}
}

}  // namespace cln
