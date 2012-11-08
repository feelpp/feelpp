// asinh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"
#include "cln/real.h"

namespace cln {

// Methode:
// Wert und Branch Cuts nach der Formel CLTL2, S. 313:
//   arsinh(z) = log(z+sqrt(1+z^2))
// z=x+iy, Ergebnis u+iv.
// Falls x=0 und y=0: u=0, v=0.
// Falls x=0: arsinh(iy) = i arcsin(y).
//   y rational ->
//     Bei y=1: u = 0, v = pi/2.
//     Bei y=1/2: u = 0, v = pi/6.
//     Bei y=0: u = 0, v = 0.
//     Bei y=-1/2: u = 0, v = -pi/6.
//     Bei y=-1: u = 0, v = -pi/2.
//     Sonst y in Float umwandeln.
//   e := Exponent aus (decode-float y), d := (float-digits y)
//   Bei y=0.0 oder e<=-d/2 liefere u = 0, v = y
//     (denn bei e<=-d/2 ist y^2/3 < y^2/2 < 2^(-d)/2 = 2^(-d-1), also
//     1 <= asin(y)/y < 1+y^2/3 < 1+2^(-d-1) < 1+2^(-d),
//     also ist asin(y)/y, auf d Bits gerundet, gleich 1.0).
//   Berechne 1-y^2.
//   Bei y>1 liefere  u = ln(y+sqrt(y^2-1)), v = pi/2.
//   Bei y<-1 liefere  u = -ln(|y|+sqrt(|y|^2-1)), v = -pi/2.
//   Bei |y|<=1 liefere  u = 0, v = atan(X=sqrt(1-y^2),Y=y).
// Falls y=0:
//   x rational -> x in Float umwandeln.
//   |x|<1/2: u = atanh(x/sqrt(1+x^2)),
//   x>=1/2: u = ln(x+sqrt(1+x^2)),
//   x<=-1/2: u = -ln(-x+sqrt(1+x^2)).
//   v = 0.
// Sonst:
//   z in Bild(sqrt) -> log(sqrt(1+z^2)+z) = (!) = 2 artanh(z/(1+sqrt(1+z^2))).
//   z nicht in Bild(sqrt) ->
//     arsinh(z) = -arsinh(-z).
//     (Denn arsinh(z)+arsinh(-z) == log((z+sqrt(1+z^2))(-z+sqrt(1+z^2)))
//           = log((1+z^2)-z^2) = log(1) = 0 mod 2 pi i, und links ist
//      der Imaginärteil betragsmäßig <=pi.)
//     Also arsinh(z) = -arsinh(-z) = - 2 artanh(-z/(1+sqrt(1+z^2)))
//          = (wegen -artanh(-w) = artanh(w)) = 2 artanh(z/(1+sqrt(1+z^2))).
// Real- und Imaginärteil des Ergebnisses sind Floats, außer wenn z reell oder
// rein imaginär ist.

// Um für zwei Zahlen u,v mit u^2-v^2=1 und u,v beide in Bild(sqrt)
// (d.h. Realteil>0.0 oder Realteil=0.0 und Imaginärteil>=0.0)
// log(u+v) zu berechnen:
//               log(u+v) = 2 artanh(v/(u+1))                            (!)
// (Beweis: 2 artanh(v/(u+1)) = log(1+(v/(u+1))) - log(1-(v/(u+1)))
//  = log((1+u+v)/(u+1)) - log((1+u-v)/(u+1)) == log((1+u+v)/(1+u-v))
//  = log(u+v) mod 2 pi i, und beider Imaginärteil ist > -pi und <= pi.)

inline const cl_C_R _asinh (const cl_N& z)
{
	if (realp(z)) {
		DeclareType(cl_R,z);
		return asinh(z,0);
	} else {
		DeclareType(cl_C,z);
		return asinh(realpart(z),imagpart(z));
	}
}

const cl_N asinh (const cl_N& z)
{
	var cl_C_R u_v = _asinh(z);
	var cl_R& u = u_v.realpart;
	var cl_R& v = u_v.imagpart;
	return complex(u,v);
}

}  // namespace cln
