// atanh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

// Methode:
// Wert und Branch Cuts nach der Formel CLTL2, S. 315:
//   artanh(z) = (log(1+z)-log(1-z)) / 2
// Sei z=x+iy, Ergebnis u+iv.
// Falls x=0 und y=0: u=0, v=0.
// Falls x=0: u = 0, v = atan(X=1,Y=y).
// Falls y=0:
//   x rational -> x in Float umwandeln.
//   |x|<1/2: u = atanh(x), v = 0.
//   |x|>=1/2: (1+x)/(1-x) errechnen,
//             =0 -> Error,
//             >0 (also |x|<1) -> u = 1/2 log((1+x)/(1-x)), v = 0.
//             <0 (also |x|>1) -> u = 1/2 log(-(1+x)/(1-x)),
//                                v = (-pi/2 für x>1, pi/2 für x<-1).
// Sonst:
//   1+x und 1-x errechnen.
//   x und y in Floats umwandeln.
//   |4x| und 1+x^2+y^2 errechnen,
//   |4x| < 1+x^2+y^2 -> u = 1/2 atanh(2x/(1+x^2+y^2)),
//   |4x| >= 1+x^2+y^2 -> u = 1/4 ln ((1+x^2+y^2)+2x)/((1+x^2+y^2)-2x)
//                        oder besser (an der Singularität: |x|-1,|y| klein):
//                        u = 1/4 ln ((1+x)^2+y^2)/((1-x)^2+y^2).
//   v = 1/2 atan(X=(1-x)(1+x)-y^2,Y=2y) * (-1 falls Y=0.0 und X<0.0 und x>=0.0,
//                                          1 sonst)
// Ergebnis ist reell nur, wenn z reell.
// Real- und Imaginärteil des Ergebnisses sind Floats, außer wenn z reell oder
// rein imaginär ist.

inline const cl_C_R _atanh (const cl_N& z)
{
	if (realp(z)) {
		DeclareType(cl_R,z);
		return atanh(z,0);
	} else {
		DeclareType(cl_C,z);
		return atanh(realpart(z),imagpart(z));
	}
}

const cl_N atanh (const cl_N& z)
{
	var cl_C_R u_v = _atanh(z);
	var cl_R& u = u_v.realpart;
	var cl_R& v = u_v.imagpart;
	return complex(u,v);
}

}  // namespace cln
