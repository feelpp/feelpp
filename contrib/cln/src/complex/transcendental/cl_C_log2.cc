// log().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"
#include "real/cl_R.h"
#include "base/cl_N.h"

namespace cln {

const cl_N log (const cl_N& a, const cl_N& b)
{
// Methode:
// (log a b) =
//   falls b reell, >0:
//     (complex (/ (log (abs a)) (log b)) (/ (phase a) (log b))), genauer:
//     falls a reell, >0: bekannt
//     falls (= a 0): Error
//     sonst: (phase a) errechnen, ein Float.
//            b (falls rational) ins selbe Float-Format umwandeln,
//            Imaginärteil := (/ (phase a) (log dieses_b)).
//            Falls a rational: (log (abs a) b).
//            Falls a komplex mit rationalem Real- und Imaginärteil,
//              Betragsquadrat  (expt (abs a) 2)  exakt ausrechnen als
//              (+ (expt (realpart a) 2) (expt (imagpart a) 2)).
//              Setze  Realteil := (/ (log Betragsquadrat b) 2).
//              [Eventuell wird hierbei (log b) ein zweites Mal ausgerechnet,
//               aber dies sowieso nur in Single-Precision.]
//            Sonst bilde (abs a), ein Float, und (log (abs a)), ein Float,
//              wandle b (falls rational) ins selbe Float-Format um,
//              setze  Realteil := (/ (log (abs a)) (log dieses_b)).
//   sonst: (/ (log a) (log b))
	if (realp(b)) {
	    DeclareType(cl_R,b);
	    if (plusp(b)) {
		// b ist reell und >0
		if (realp(a)) {
			DeclareType(cl_R,a);
			if (plusp(a))
				// a und b sind beide reell und >0
				return log(a,b);
		}
		// b ist reell und >0, a aber nicht.

		// Imaginärteil (/ (phase a) (log b)) errechnen:
		var cl_F im;
		{
			var cl_R angle = phase(a);
			if (eq(angle,0)) // = Fixnum 0 <==> (= a 0) -> Error
				{ throw division_by_0_exception(); }
		 {	DeclareType(cl_F,angle);
			var cl_F bf = cl_somefloat(b,angle); // (float b)
			im = angle / ln(bf);
		}}

		// Realteil (/ (log (abs a)) (log b)) errechnen:
		var cl_R re;
		if (realp(a)) {
			DeclareType(cl_R,a);
			if (rationalp(a)) {
				// a rational -> (log (abs a) b) errechnen:
				re = log(abs(a),b); // NB: (abs a) > 0
				goto re_ok;
			}
		} else {
			DeclareType(cl_C,a);
			if (rationalp(realpart(a)) && rationalp(imagpart(a))) {
				// a komplex mit rationalem Real- und Imaginärteil a1,a2
				var const cl_R& a1 = realpart(a);
				var const cl_R& a2 = imagpart(a);
				re = log(square(a1)+square(a2),b) / 2;
				goto re_ok;
			}
		}
		// Keine Chance für rationalen Realteil.
		{
			var cl_F abs_a = The(cl_F)(abs(a));
			var cl_F log_abs_a = ln(abs_a);
			var cl_F bf = cl_somefloat(b,log_abs_a); // (float b)
			re = log_abs_a / ln(bf);
		}
		re_ok:

		return complex_C(re,im);
	    }
	}

	// normaler komplexer Fall
	return log(a) / log(b);
}

}  // namespace cln
