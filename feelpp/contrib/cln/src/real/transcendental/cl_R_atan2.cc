// atan().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "cln/float.h"
#include "float/transcendental/cl_F_tran.h"
#include "base/cl_N.h"
#include "real/cl_R.h"

namespace cln {

const cl_R atan (const cl_R& x, const cl_R& y)
{
// Methode:
// y=0 -> bei x>0: 0 als Ergebnis,
//        bei x<0: pi als Ergebnis.
//        bei x=0: Error.
// x=0 -> bei y>0: pi/2 als Ergebnis.
//        bei y<0: -pi/2 als Ergebnis.
//        bei y=0: Error.
// Falls x und y beide rational: beide in Floats umwandeln.
// 0 <= |y| <= x  ->  atan(y/x)
// 0 <= |x| <= y  ->  pi/2 - atan(x/y)
// 0 <= |x| <= -y  ->  -pi/2 - atan(x/y)
// 0 <= |y| <= -x  ->  für y>=0: pi + atan(y/x), für y<0: -pi + atan(y/x)

	if (eq(y,0)) {
		// y=0 (exakt)
		if (zerop(x)) // x=0 -> Error
			{ throw division_by_0_exception(); }
		if (minusp(x)) // x<0 -> pi in Default-Float-Genauigkeit
			{ return pi(); }
		return 0; // x>0 -> 0
	}
	elif (eq(x,0)) {
		// x=0 (exakt)
		if (zerop(y)) // y=0 -> Error
			{ throw division_by_0_exception(); }
		if (minusp(y)) // y<0 -> -pi/2
			{ return - scale_float(pi(),-1); }
		return scale_float(pi(),-1); // y>0 -> pi/2
	} else {
		Mutable(cl_R,x);
		Mutable(cl_R,y);
		// Check special case of rational numbers:
		if (rationalp(x))
			if (rationalp(y)) {
				// x,y in Floats umwandeln:
				x = cl_float(The(cl_RA)(x));
				y = cl_float(The(cl_RA)(y));
			}
		// x,y nicht exakt =0, x/y und y/x werden Floats sein.
		if (abs(x) >= abs(y)) {
			// |x| >= |y|
			var cl_F z = atanx(The(cl_F)(y/x));
			// Division war erfolgreich, also x/=0.
			if (minusp(x))
				// x<0 -> pi bzw. -pi addieren:
				if (!minusp(y))
					// y>=0 -> atan(y/x) + pi
					return z + pi(z);
				else
					// y<0 -> atan(y/x) - pi
					return z - pi(z);
			else
				return z;
		} else {
			// |x| < |y|
			var cl_F z = atanx(The(cl_F)(x/y));
			// von pi/2 bzw. -pi/2 subtrahieren:
			if (!minusp(y))
				// y>=0 -> pi/2 - atan(x/y)
				return scale_float(pi(z),-1) - z;
			else
				// y<0 -> -pi/2 - atan(x/y)
				return - scale_float(pi(z),-1) - z;
		}
	}
}

}  // namespace cln
