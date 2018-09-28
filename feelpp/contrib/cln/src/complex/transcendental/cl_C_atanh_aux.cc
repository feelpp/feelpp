// atanh().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "complex/cl_C.h"


// Implementation.

#include "base/cl_N.h"
#include "cln/real.h"
#include "float/transcendental/cl_F_tran.h"
#include "real/cl_R.h"

/* Use the inline version of cl_float */
#include "base/cl_inline.h"
#include "real/conv/cl_F_from_R_def.cc"

namespace cln {

// Hilfsfunktion für atanh und atan: u+iv := artanh(x+iy). Liefert cl_C_R(u,v).

const cl_C_R CL_FLATTEN atanh (const cl_R& x, const cl_R& y)
{
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

	if (eq(x,0))
		// x=0 -> u=0, v=atan(X=1,Y=y) (Fall y=0 ist inbegriffen)
		return cl_C_R(0, atan(1,y));
	if (eq(y,0)) {
		var cl_F xf = cl_float_inline(x); // (float x)
		var cl_F& x = xf;
		// x Float
		if (zerop(x))
			// x=0.0 -> x als Ergebnis
			return cl_C_R(x, 0);
		if (float_exponent(x) < 0)
			// Exponent e<0, also |x|<1/2
			return cl_C_R(atanhx(x), 0);
		// e>=0, also |x|>=1/2
		var cl_F xx_den = 1 - x;
		var cl_F xx = (1 + x) / xx_den; // (1+x)/(1-x)
		var cl_R v;
		if (!minusp(xx)) {
			if (zerop(xx))
				{ throw division_by_0_exception(); }
			v = 0;
		} else {
			// (1+x)/(1-x) < 0 -> Betrag nehmen, Imaginärteil berechnen:
			xx = - xx;
			v = scale_float(pi(),-1); // (scale-float pi -1) = pi/2
			if (minusp(xx_den))
				// 1-x<0 -> dann -pi/2
				v = -v;
		}
		// ln bilden, durch 2
		return cl_C_R(scale_float(ln(xx),-1), v);
	}
	var cl_R _1_plus_x = 1+x;
	var cl_R _1_minus_x = 1-x;
	// x und y in Floats umwandeln: (Diese Fallunterscheidung ist
	// symmetrisch in x und y, auch wenn's nicht so aussieht.)
	var cl_F xf;
	var cl_F yf;
	if (rationalp(x)) {
		DeclareType(cl_RA,x);
		yf = cl_float_inline(y);
		xf = cl_float(x,yf);
	} else {
		DeclareType(cl_F,x);
		xf = x;
		yf = cl_somefloat(y,xf);
	}
	var cl_F yf_2 = square(yf);
	var cl_F u;
	{
		var cl_F temp1 = abs(scale_float(xf,2)); // |4x|
		var cl_F temp2 = 1 + (square(xf) + yf_2); // 1+x^2+y^2
		if (temp1 < temp2) // |4x| < 1+x^2+y^2 ?
			// u = 1/2 atanh(2x/(1+x^2+y^2))
			u = scale_float(atanhx(scale_float(xf,1)/temp2),-1);
		else {
			// u = 1/4 ln ((1+x)^2+y^2)/((1-x)^2+y^2)
			var cl_F num = _1_plus_x*_1_plus_x + yf_2; // (1+x)^2+y^2, ein Float >=0
			var cl_F den = _1_minus_x*_1_minus_x + yf_2; // (1-x)^2+y^2, ein Float >=0
			if (zerop(den))
				{ throw division_by_0_exception(); }
			u = scale_float(ln(num/den),-2);
		}
	}
	var cl_F v;
	{
		var cl_F X = _1_plus_x*_1_minus_x-yf_2;
		var cl_F Y = scale_float(yf,1);
		v = atan(X,Y); // atan(X=(1-x)(1+x)-y^2,Y=2y), ein Float
		if (minusp(X) && !minusp(x) && zerop(Y))
			v = -v;
		v = scale_float(v,-1); // 1/2 * atan(...) * +-1
	}
	return cl_C_R(u,v);
}

}  // namespace cln
