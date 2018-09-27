// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

cl_signean compare (const cl_RA& r, const cl_RA& s)
{
// Methode:
// r,s Integer -> klar
// r<0, s>=0 -> r<s.
// r>=0, s<0 -> r>s.
// r Integer, s Ratio: r=a, s=b/c. Vergleiche a*c und b.
// r Ratio, s Integer: r=a/b, s=c. Vergleiche a und b*c.
// r,s Ratios: r=a/b, s=c/d. Vergleiche a*d und b*c.
	// 1. Schritt: Test, ob beides Integers:
	if (integerp(r) && integerp(s)) {
		DeclareType(cl_I,r);
		DeclareType(cl_I,s);
		return compare(r,s);
	}
	// r,s nicht beide Integers.
	// 2. Schritt: Test, ob die Vorzeichen bereits das Ergebnis hergeben:
	if (minusp(r)) {
		if (!minusp(s))
			return signean_minus; // r<0, s>=0 -> r<s
	} else {
		if (minusp(s))
			return signean_plus; // r>=0, s<0 -> r>s
	}
	// r,s haben gleiches Vorzeichen.
	// 3. Schritt: Fallunterscheidung nach Typen
	if (integerp(r)) {
		DeclareType(cl_I,r);
		DeclareType(cl_RT,s);
		// r Integer, s Ratio: r=a, s=b/c. Vergleiche a*c und b.
		var const cl_I& a = r;
		var const cl_I& b = numerator(s);
		var const cl_I& c = denominator(s);
		return compare(a*c,b);
	}
	elif (integerp(s)) {
		DeclareType(cl_I,s);
		DeclareType(cl_RT,r);
		// r Ratio, s Integer: r=a/b, s=c. Vergleiche a und b*c.
		var const cl_I& a = numerator(r);
		var const cl_I& b = denominator(r);
		var const cl_I& c = s;
		return compare(a,b*c);
	}
	else {
		DeclareType(cl_RT,r);
		DeclareType(cl_RT,s);
		// r,s Ratios: r=a/b, s=c/d. Vergleiche a*d und b*c.
		var const cl_I& a = numerator(r);
		var const cl_I& b = denominator(r);
		var const cl_I& c = numerator(s);
		var const cl_I& d = denominator(s);
		return compare(a*d,b*c);
	}
}
// Beschleunigung durch Konversion zu Short-Floats diese zuerst vergleichen??

}  // namespace cln
