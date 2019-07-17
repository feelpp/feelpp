// acos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"
#include "real/cl_R.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"
#include "cln/float.h"

namespace cln {

inline const cl_F pi (const cl_R& v)
{
	if (rationalp(v))
		return pi();
	else {
		DeclareType(cl_F,v);
		return pi(v);
	}
}

const cl_N acos (const cl_N& z)
{
// Methode:
// Wert und Branch Cuts nach der Formel CLTL2, S. 312:
//   arccos(z) = log(z+i*sqrt(1-z^2))/i = pi/2 - arcsin(z)
// Sei z=x+iy.
// Falls y=0:
//   Falls x rational:
//     Bei x=1: Ergebnis 0.
//     Bei x=1/2: Ergebnis pi/3.
//     Bei x=0: Ergebnis pi/2.
//     Bei x=-1/2: Ergebnis 2pi/3.
//     Bei x=-1: Ergebnis pi.
//     Sonst x in Float umwandeln.
//   Falls x>1: Ergebnis i ln(x+sqrt(x^2-1)).
// Sonst errechne u+iv = arsinh(-y+ix) wie oben, Ergebnis (pi/2-v)+iu.

	cl_C_R u_v;
	if (realp(z)) {
		DeclareType(cl_R,z);
		// y=0
		var const cl_R& x = z;
		var cl_F xf;
		if (rationalp(x)) {
			DeclareType(cl_RA,x);
			// x rational
			if (integerp(x)) {
				DeclareType(cl_I,x);
				// x Integer
				if (eq(x,0)) // x=0 -> Ergebnis pi/2
					return scale_float(pi(),-1);
				if (eq(x,1)) // x=1 -> Ergebnis 0
					return 0;
				if (eq(x,-1)) // x=-1 -> Ergebnis pi
					return pi();
				xf = cl_float(x);
			} else {
				DeclareType(cl_RT,x);
				// x Ratio
				if (eq(denominator(x),2)) { // Nenner = 2 ?
					if (eq(numerator(x),1)) // x=1/2 -> Ergebnis pi/3
						return pi()/3;
					if (eq(numerator(x),-1)) // x=-1/2 -> Ergebnis 2pi/3
						return scale_float(pi(),1)/3;
				}
				xf = cl_float(x);
			}
		} else {
			DeclareType(cl_F,x);
			xf = x;
		}
		// x Float
	 {	var cl_F& x = xf;
		if (cl_I(1) < x)
			// x>1
			return complex_C(0,ln(x+sqrt(square(x)-1)));
		u_v = asinh(0,x);
	}} else {
		DeclareType(cl_C,z);
		u_v = asinh(-imagpart(z),realpart(z));
	}
	var cl_R& u = u_v.realpart;
	var cl_R& v = u_v.imagpart;
	var cl_F archimedes = pi(v); // pi im Float-Format von v
	return complex(scale_float(archimedes,-1)-v,u); // (pi/2-v)+iu
}

}  // namespace cln
