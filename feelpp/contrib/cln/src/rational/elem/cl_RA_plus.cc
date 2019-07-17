// binary operator +

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_RA operator+ (const cl_RA& r, const cl_RA& s)
{
// Methode (vgl. [Buchberger, Collins, Loos: Computer Algebra, S.200-201])
// r,s beide Integers -> klar.
// r=a/b, s=c -> Ergebnis (a+b*c)/b
//   (mit b>1 und ggT(a+b*c,b) = ggT(a,b) = 1)
//   Bei c=0 direkt r als Ergebnis.
// r=a, s=c/d -> Ergebnis (a*d+c)/d
//   (mit d>1 und ggT(a*d+c,d) = ggT(c,d) = 1)
//   Bei a=0 direkt s als Ergebnis.
// r=a/b, s=c/d:
//   g:=ggT(b,d)>0.
//   Falls g=1:
//     Ergebnis (a*d+b*c)/(b*d),
//     (mit b*d>1 wegen b>1, d>1, und
//      ggT(a*d+b*c,b*d) = 1
//      wegen ggT(a*d+b*c,b) = ggT(a*d,b) = 1 (wegen ggT(a,b)=1 und ggT(d,b)=1)
//      und   ggT(a*d+b*c,d) = ggT(b*c,d) = 1 (wegen ggT(b,d)=1 und ggT(c,d)=1)
//     )
//   Sonst b' := b/g, d' := d/g. e := a*d'+b'*c, f:= b'*d = b*d'.
//   Es ist g = ggT(g*b',g*d') = g*ggT(b',d'), also ggT(b',d')=1.
//   Es ist r+s = (a*d+b*c)/(b*d) = (nach Kürzen mit g) e/f.
//   Außerdem:
//     ggT(a,b') teilt ggT(a,b)=1, also ggT(a,b')=1. Mit ggT(d',b')=1 folgt
//     1 = ggT(a*d',b') = ggT(a*d'+b'*c,b') = ggT(e,b').
//     ggT(c,d') teilt ggT(c,d)=1, also ggT(c,d')=1. Mit ggT(b',d')=1 folgt
//     1 = ggT(b'*c,d') = ggT(a*d'+b'*c,d') = ggT(e,d').
//     Daher ist ggT(e,f) = ggT(e,b'*d'*g) = ggT(e,g).
//   Errechne daher h=ggT(e,g).
//   Bei h=1 ist e/f das Ergebnis (mit f>1, da d>1, und ggT(e,f)=1),
//   sonst ist (e/h)/(f/h) das Ergebnis.
	if (integerp(s)) {
		// s ist Integer
		DeclareType(cl_I,s);
		if (eq(s,0)) { return r; } // s=0 -> r als Ergebnis
		if (integerp(r)) {
			// beides Integers
			DeclareType(cl_I,r);
			return r+s;
		} else {
			DeclareType(cl_RT,r);
			var const cl_I& a = numerator(r);
			var const cl_I& b = denominator(r);
			var const cl_I& c = s;
			// r = a/b, s = c.
			return I_I_to_RT(a+b*c,b);
		}
	} else {
		// s ist Ratio
		DeclareType(cl_RT,s);
		if (integerp(r)) {
			// r ist Integer
			DeclareType(cl_I,r);
			if (eq(r,0)) { return s; } // r=0 -> s als Ergebnis
			var const cl_I& a = r;
			var const cl_I& c = numerator(s);
			var const cl_I& d = denominator(s);
			// r = a, s = c/d.
			return I_I_to_RT(a*d+c,d);
		} else {
			// r,s beide Ratios
			DeclareType(cl_RT,r);
			var const cl_I& a = numerator(r);
			var const cl_I& b = denominator(r);
			var const cl_I& c = numerator(s);
			var const cl_I& d = denominator(s);
			var cl_I g = gcd(b,d); // g = ggT(b,d) >0 bilden
			if (eq(g,1))
				// g=1 -> Ergebnis (a*d+b*c)/(b*d)
				return I_I_to_RT(a*d+b*c,b*d);
			// g>1
			var cl_I bp = exquopos(b,g); // b' := b/g (b,g>0)
			var cl_I dp = exquopos(d,g); // d' := d/g (d,g>0)
			var cl_I e = a*dp+bp*c; // e := a*d'+b'*c
			var cl_I f = bp*d; // f := b'*d
			var cl_I h = gcd(e,g); // h := ggT(e,g)
			if (eq(h,1))
				// h=1
				return I_I_to_RT(e,f);
			// h>1
			return I_I_to_RA(exquo(e,h),exquopos(f,h)); // (e/h)/(f/h) als Ergebnis
		}
	}
}

}  // namespace cln
