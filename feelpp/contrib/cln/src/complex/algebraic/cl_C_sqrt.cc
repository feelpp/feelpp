// sqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

const cl_N sqrt (const cl_N& x)
{
// Methode:
// x reell -> Für x>=0 klar, für x<0: sqrt(-x)*i.
// x=a+bi ->
//   Bestimme r=abs(x)=sqrt(a*a+b*b).
//   Falls a>=0: Setze c:=sqrt((r+a)/2), d:=(b/(2*c) falls c>0, c falls c=0).
//   Falls a<0: Setze d:=sqrt((r-a)/2)*(1 falls b>=0, -1 falls b<0), c:=b/(2*d).
//   Damit ist c>=0, 2*c*d=b, c*c=(r+a)/2, d*d=(r-a)/2, c*c-d*d=a, c*c+d*d=r,
//   also c+di die gesuchte Wurzel.
	if (realp(x)) {
		DeclareType(cl_R,x);
		if (!minusp(x))
			return sqrt(x);
		else
			return complex_C(0,sqrt(-x));
	} else {
		DeclareType(cl_C,x);
		var const cl_R& a = realpart(x);
		var const cl_R& b = imagpart(x);
		var cl_R r = cl_hypot(a,b); // r = (abs x)
		if (!minusp(a)) {
			// a>=0
			var cl_R c = sqrt((r+a)/2);
			var cl_R d = (!zerop(c) ? b/(2*c) : c);
			return complex_C(c,d);
		} else {
			var cl_R d = sqrt((r-a)/2);
			if (minusp(b))
				d = -d;
			var cl_R c = b/(2*d);
			return complex_C(c,d);
		}
	}
}

}  // namespace cln
