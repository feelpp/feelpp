// binary operator /

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N operator/ (const cl_N& x, const cl_N& y)
{
// Methode:
// x,y beide reell -> klar.
// x=a+bi, y=c reell -> (a/c)+(b/c)i
// y komplex -> (* x (/ y))
	if (realp(y)) {
		DeclareType(cl_R,y);
		if (realp(x)) {
			DeclareType(cl_R,x);
			// x,y beide reell
			return x/y;
		} else {
			DeclareType(cl_C,x);
			// x komplex: x=a+bi, y=c
			var const cl_R& a = realpart(x);
			var const cl_R& b = imagpart(x);
			var const cl_R& c = y;
			return complex(a/c,b/c);
		}
	} else {
		DeclareType(cl_C,y);
		// y komplex
		return x * recip(y);
	}
}

}  // namespace cln
