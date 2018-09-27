// asinh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

#include "cln/real.h"
#include "float/transcendental/cl_F_tran.h"
#include "real/cl_R.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"
#include "cln/float.h"

/* Use the inline version of cl_float */
#include "base/cl_inline.h"
#include "real/conv/cl_F_from_R_def.cc"

namespace cln {

// Hilfsfunktion für asinh und asin: u+iv := arsinh(x+iy). Liefert cl_C_R(u,v).

const cl_C_R asinh (const cl_R& x, const cl_R& y)
{
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

	if (eq(x,0)) {
		// x=0
		var cl_F yf;
		if (rationalp(y)) {
			DeclareType(cl_RA,y);
			// y rational
			if (eq(y,0)) // x=0, y=0 -> u=0, v=0
				return cl_C_R(0,0);
			if (integerp(y)) {
				DeclareType(cl_I,y);
				// y Integer
				if (eq(y,1)) // x=0, y=1 -> v = pi/2
					return cl_C_R(0,scale_float(pi(),-1));
				if (eq(y,-1)) // x=0, y=-1 -> v = -pi/2
					return cl_C_R(0,-scale_float(pi(),-1));
				yf = cl_float_inline(y); // y in Float umwandeln
			} else {
				DeclareType(cl_RT,y);
				// y Ratio
				if (eq(denominator(y),2)) { // Nenner = 2 ?
					if (eq(numerator(y),1)) // x=0, y=1/2 -> v = pi/6
						return cl_C_R(0,pi()/6);
					if (eq(numerator(y),-1)) // x=0, y=-1/2 -> v = -pi/6
						return cl_C_R(0,-(pi()/6));
				}
				yf = cl_float_inline(y); // y in Float umwandeln
			}
		} else {
			DeclareType(cl_F,y);
			yf = y;
		}
		// y Float
		var cl_F& y = yf;
		if (zerop(y)) // y=0.0 -> arcsin(y) = y als Ergebnis
			return cl_C_R(0,y);
		if (float_exponent(y) <= (-(sintC)float_digits(y))>>1)
			// e <= -d/2 <==> e <= -ceiling(d/2)
			return cl_C_R(0,y);
		var cl_F temp = 1-square(y);
		if (!minusp(temp))
			// 1-y*y>=0, also |y|<=1
			// v = atan(X=sqrt(1-y*y),Y=y)
			return cl_C_R(0,atan(sqrt(temp),y));
		else {
			// 1-y*y<0, also |y|>1
			temp = sqrt(-temp); // sqrt(y*y-1)
			if (minusp(y))
				temp = temp - y;
			else
				temp = temp + y;
			// temp = sqrt(y^2-1)+|y|, ein Float >1
			var cl_F u = ln(temp); // ln(|y|+sqrt(y^2-1)), ein Float >0
			var cl_F v = scale_float(pi(),-1); // (scale-float pi -1) = pi/2
			if (!minusp(y))
				return cl_C_R(u,v); // y>1 -> v = pi/2
			else
				return cl_C_R(-u,-v); // y<-1 -> v = -pi/2, u = -ln(...)
		}
	}
	if (eq(y,0)) {
		// y=0
		var cl_F xf = cl_float_inline(x); // x in Float umwandeln
		var cl_F& x = xf;
		// x Float
		if (zerop(x))
			return cl_C_R(x,0); // x=0.0 -> u=x, v=0.
		var cl_F temp = sqrt(1+square(x)); // sqrt(1+x^2)
		if (float_exponent(x) < 0) // Exponent e (von x/=0) <0 ?
			// |x|<1/2
			return cl_C_R(atanhx(x/temp),0);
		else
			// |x|>=1/2
			if (!minusp(x))
				// x>=1
				return cl_C_R(ln(temp+x),0); // u = ln(x+sqrt(1+x^2))
			else
				// x<=-1
				return cl_C_R(-ln(temp-x),0); // u = -ln(-x+sqrt(1+x^2))
	}
	var cl_N z = complex_C(x,y); // z=x+iy
	var cl_N w = z/(1+sqrt(1+square(z))); // z/(1+sqrt(1+z^2))
	// Da z=x+iy weder reell noch rein imaginär ist, ist auch
	// w := z/(1+sqrt(1+z^2)) weder reell noch rein imaginär.
	// (Beweis: Sollte sqrt(1+z^2) rationalen Real- und Imaginärteil haben,
	// so auch z, also auch w, und die Formel z = 2w/(1-w^2) zeigt, daß dann
	// z reell oder rein imaginär sein müßte. Also hat sqrt(1+z^2) ein
	// Float als Real- oder Imaginärteil, das Betragsquadrat des Nenners
	// ist also ein Float, und da Real- und Imaginärteil von z /=0 sind,
	// sind Real- und Imaginärteil von w Floats.)
	// Daher hat dann atanh(...) Floats als Realteil u und Imaginärteil v.
 {	DeclareType(cl_C,w);
	cl_C_R u_v = atanh(realpart(w),imagpart(w));
	var cl_R& u = u_v.realpart;
	var cl_R& v = u_v.imagpart;
  {	DeclareType(cl_F,u);
	DeclareType(cl_F,v);
	return cl_C_R(scale_float(u,1),scale_float(v,1)); // u:=2*u, v:=2*v
}}}

}  // namespace cln
