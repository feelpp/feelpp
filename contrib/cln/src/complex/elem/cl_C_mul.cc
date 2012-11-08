// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N operator* (const cl_N& x, const cl_N& y)
{
// Methode:
// x,y beide reell -> klar.
// x=a, y=b+ci -> (a*b)+(a*c)i
// x=a+bi, y=c -> (a*c)+(b*c)i
// x=a+bi, y=c+di -> (a*c-b*d)+(a*d+b*c)i
	if (realp(x)) {
		DeclareType(cl_R,x);
		if (realp(y)) {
			DeclareType(cl_R,y);
			return x*y;
		} else {
			DeclareType(cl_C,y);
			// x=a, y=b+ci
			var const cl_R& a = x;
			var const cl_R& b = realpart(y);
			var const cl_R& c = imagpart(y);
			return complex(a*b,a*c);
		}
	} else {
		DeclareType(cl_C,x);
		if (realp(y)) {
			DeclareType(cl_R,y);
			// x=a+bi, y=c
			var const cl_R& a = realpart(x);
			var const cl_R& b = imagpart(x);
			var const cl_R& c = y;
			return complex(a*c,b*c);
		} else {
			DeclareType(cl_C,y);
			// x=a+bi, y=c+di
			var const cl_R& a = realpart(x);
			var const cl_R& b = imagpart(x);
			var const cl_R& c = realpart(y);
			var const cl_R& d = imagpart(y);
			return complex(a*c-b*d,a*d+b*c);
		}
	}
}

}  // namespace cln
