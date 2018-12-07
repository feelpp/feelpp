// expt().

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
#include "integer/cl_I.h"
#include "base/cl_N.h"

namespace cln {

// Methode:
// Falls y rational:
//   Falls y Integer:
//     Falls y=0: Ergebnis 1,
//       [Nach CLTL folgendermaßen:
//         x reell:
//           x rational -> Fixnum 1
//           x Float -> (float 1 x)
//         x komplex:
//           x komplex rational -> Fixnum 1
//           sonst: #C(1.0 0.0) im Float-Format des Real- bzw. Imaginärteils von x
//       ]
//     Falls x rational oder komplex rational oder |y| klein:
//       x^|y| durch wiederholtes Quadrieren und Multiplizieren und evtl.
//       Kehrwert-Bilden ermitteln.
//     Sonst wie bei 'y Float'.
//   Falls y Ratio m/n:
//     Es gilt (expt x m/n) = (expt (expt x 1/n) m).
//     Falls x in Q(i) liegt (also rational oder komplex rational ist):
//       Sollte x^(m/n) in Q(i) liegen, so auch eine n-te Wurzel x^(1/n)
//       (und bei n=2 oder n=4 damit auch alle n-ten Wurzeln x^(1/n) ).
//       Falls x rational >=0: n-te Wurzel aus x nehmen. Ist sie rational,
//         deren m-te Potenz als Ergebnis.
//       Falls x rational <=0 oder komplex rational und n Zweierpotenz:
//         n-te Wurzel aus x nehmen (mehrfaches sqrt). Ist sie rational oder
//         komplex rational, deren m-te Potenz als Ergebnis.
//         [Beliebige n betrachten!??]
//     Falls n Zweierpotenz und |m|,n klein: n-te Wurzel aus x nehmen
//       (mehrfaches sqrt), davon die m-te Potenz durch wiederholtes
//       Quadrieren und Multiplizieren und evtl. Kehrwert-Bilden.
//     Sonst wie bei 'y Float'.
// Falls y Float oder komplex:
//   Falls (zerop x):
//     Falls Realteil von y >0 :
//       liefere 0.0 falls x und y reell, #C(0.0 0.0) sonst.
//     Sonst Error.
//   Falls y=0.0:
//     liefere 1.0 falls x und y reell, #C(1.0 0.0) sonst.
//   Sonst: (exp (* (log x) y))
// Das Ergebnis liegt in Q(i), falls x in Q(i) liegt und 4y ein Integer ist.??
// Genauigkeit erhöhen, log2(|y|) Bits mehr??
// Bei x oder y rational und der andere Long-Float: bitte kein Single-Float!??

// Liefert x^0.
inline const cl_N expt_0 (const cl_N& x)
{
#ifdef STRICT_CLTL
	// y=0 -> 1 im Format von x.
	if (realp(x)) {
		DeclareType(cl_R,x);
		if (rationalp(x)) {
			DeclareType(cl_RA,x);
			return 1;
		} else {
			DeclareType(cl_F,x);
			return cl_float(1,x);
		}
	} else {
		DeclareType(cl_C,x);
		var cl_R f = contagion(realpart(x),imagpart(x));
		if (rationalp(f)) {
			DeclareType(cl_RA,f);
			return 1;
		} else {
			DeclareType(cl_F,f);
			// #C(1.0 0.0)
			return complex_C(cl_float(1,f),cl_float(0,f));
		}
	}
#else
	// Exponent exakt 0 -> Ergebnis exakt 1
	unused x;
	return 1;
#endif
}

inline const cl_R contagion (const cl_N& x)
{
	if (realp(x)) {
		DeclareType(cl_R,x);
		return x;
	} else {
		DeclareType(cl_C,x);
		return contagion(realpart(x),imagpart(x));
	}
}

const cl_N expt (const cl_N& x, const cl_N& y)
{
	if (realp(y)) {
	    DeclareType(cl_R,y);
	    if (rationalp(y)) {
		DeclareType(cl_RA,y);
		// y rational
		if (integerp(y)) {
			DeclareType(cl_I,y);
			// y Integer
			if (eq(y,0))
				return expt_0(x); // Liefere 1
			if (fixnump(y)) // |y| klein ?
				return expt(x,y); // exakt ausrechnen
			if (realp(x)) {
				DeclareType(cl_R,x);
				if (rationalp(x)) {
					DeclareType(cl_RA,x);
					return expt(x,y); // exakt ausrechnen
				}
			} else {
				DeclareType(cl_C,x);
				if (rationalp(realpart(x)) && rationalp(imagpart(x)))
					return expt(x,y); // exakt ausrechnen
			}
			// x nicht exakt und |y| groß
		} else {
			DeclareType(cl_RT,y);
			// y Ratio
			var const cl_I& m = numerator(y);
			var const cl_I& n = denominator(y);
			if (realp(x)) {
				DeclareType(cl_R,x);
				if (rationalp(x)) {
					DeclareType(cl_RA,x);
					if (minusp(x))
						goto complex_rational;
					// x rational >=0
					var cl_RA w;
					if (rootp(x,n,&w)) // Wurzel rational?
						return expt(w,m);
				}
			} else {
				DeclareType(cl_C,x);
				if (rationalp(realpart(x)) && rationalp(imagpart(x)))
					goto complex_rational;
			}
			if (false) {
				complex_rational:
				// x in Q(i)
				var uintC k = power2p(n);
				if (k) {
					// n Zweierpotenz = 2^(k-1). n>1, also k>1
					Mutable(cl_N,x);
					k--;
					do { x = sqrt(x); }
					   while (--k > 0);
					return expt(x,m);
				}
			}
			if (fixnump(m) && fixnump(n)) { // |m| und n klein?
				var uintV _n = FN_to_UV(n);
				if ((_n & (_n-1)) == 0) { // n Zweierpotenz?
					Mutable(cl_N,x);
					until ((_n = _n >> 1) == 0)
					      { x = sqrt(x); }
					return expt(x,m);
				}
			}
		}
	    }
	}
	// allgemeiner Fall (z.B. y Float oder komplex):
	if (zerop(x)) { // x=0.0 ?
		if (zerop(y))  // y=0.0?
			return expt_0(x); // Liefere 1
		if (rationalp(realpart(y))) // Realteil von y >0 exakt.
			return 0;
		if (!plusp(realpart(y))) // Realteil von y <=0 ?
			throw division_by_0_exception();
		else {
			var cl_R f = contagion(contagion(x),contagion(y));
			// ein Float, da sonst x = Fixnum 0 gewesen wäre
		    {	DeclareType(cl_F,f);
			var cl_F f0 = cl_float(0,f);
			return complex_C(f0,f0); // #C(0.0 0.0)
		    }
		}
	}
	return exp(log(x)*y);
}

}  // namespace cln
