// equal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

bool equal (const cl_N& x, const cl_N& y)
{
// Methode:
// Falls beide reell, klar.
// Falls x reell, y komplex: (= x (realpart y)) und (zerop (imagpart y)).
// Falls x komplex, y reell: analog
// Falls beide komplex: Realteile und Imaginärteile jeweils gleich?
	if (realp(x)) {
		DeclareType(cl_R,x);
		if (realp(y)) {
			DeclareType(cl_R,y);
			// x,y beide reell
			return equal(x,y);
		} else {
			DeclareType(cl_C,y);
			// x reell, y komplex
			if (!zerop(imagpart(y)))
				return false;
			return equal(x,realpart(y));
		}
	} else {
		DeclareType(cl_C,x);
		if (realp(y)) {
			DeclareType(cl_R,y);
			// x komplex, y reell
			if (!zerop(imagpart(x)))
				return false;
			return equal(realpart(x),y);
		} else {
			DeclareType(cl_C,y);
			// x,y beide komplex
			if (!equal(realpart(x),realpart(y)))
				return false;
			if (!equal(imagpart(x),imagpart(y)))
				return false;
			return true;
		}
	}
}

}  // namespace cln
