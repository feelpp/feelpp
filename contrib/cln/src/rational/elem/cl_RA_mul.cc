// binary operator *

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"
#include "integer/cl_I.h"

namespace cln {

const cl_RA operator* (const cl_RA& r, const cl_RA& s)
{
// Methode (vgl. [Buchberger, Collins, Loos: Computer Algebra, S.201])
// r,s beide Integers -> klar.
// r=a/b, s=c ->
//   Bei c=0 Ergebnis 0.
//   g:=ggT(b,c).
//   Falls g=1: Ergebnis (a*c)/b (mit b>1, ggT(a*c,b)=1).
//   Sonst: b':=b/g, c':=c/g, Ergebnis (a*c')/b' (mit ggT(a*c',b')=1).
// r=a, s=c/d analog.
// r=a/b, s=c/d ->
//   g:=ggT(a,d), h:=ggT(b,c).
//   a':=a/g, d':=d/g (nur bei g>1 bedeutet das Rechnung).
//   b':=b/h, c':=c/h (nur bei h>1 bedeutet das Rechnung).
//   Ergebnis ist = (a'*c')/(b'*d').
	if (integerp(s)) {
		// s Integer
		DeclareType(cl_I,s);
		if (integerp(r)) {
			// beides Integer
			DeclareType(cl_I,r);
			return r*s;
		} else {
			DeclareType(cl_RT,r);
			var const cl_I& a = numerator(r);
			var const cl_I& b = denominator(r);
			var const cl_I& c = s;
			// r=a/b, s=c, bilde a/b * c.
			if (zerop(c))
				{ return 0; } // c=0 -> Ergebnis 0
			var cl_I g = gcd(b,c);
			if (eq(g,1))
				// g=1
				return I_I_to_RT(a*c,b); // (a*c)/b
			else
				// g>1
				return I_I_to_RA(a*exquo(c,g),exquopos(b,g)); // (a*(c/g))/(b/g)
		}
	} else {
		// s ist Ratio
		DeclareType(cl_RT,s);
		if (integerp(r)) {
			// r Integer
			DeclareType(cl_I,r);
			var const cl_I& a = r;
			var const cl_I& b = numerator(s);
			var const cl_I& c = denominator(s);
			// r=a, s=b/c, bilde a * b/c.
			if (zerop(a))
				{ return 0; } // a=0 -> Ergebnis 0
			var cl_I g = gcd(a,c);
			if (eq(g,1))
				// g=1
				return I_I_to_RT(a*b,c); // (a*b)/c
			else
				// g>1
				return I_I_to_RA(exquo(a,g)*b,exquopos(c,g)); // ((a/g)*b)/(c/g)
		} else {
			// r,s beide Ratios
			DeclareType(cl_RT,r);
			var const cl_I& a = numerator(r);
			var const cl_I& b = denominator(r);
			var const cl_I& c = numerator(s);
			var const cl_I& d = denominator(s);
			var cl_I ap, dp;
			{
				var cl_I g = gcd(a,d);
				if (eq(g,1))
					{ ap = a; dp = d; }
				else
					{ ap = exquo(a,g); dp = exquopos(d,g); }
			}
			var cl_I cp, bp;
			{
				var cl_I h = gcd(b,c);
				if (eq(h,1))
					{ cp = c; bp = b; }
				else
					{ cp = exquo(c,h); bp = exquopos(b,h); }
			}
			return I_I_to_RA(ap*cp,bp*dp); // (a'*c')/(b'*d')
		}
	}
}

}  // namespace cln
