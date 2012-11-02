// rational().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "cln/integer.h"
#include "rational/cl_RA.h"

namespace cln {

const cl_RA rational (const cl_F& x)
{
  // Methode:
  // Der mathematische Wert eines Float ist, wenn INTEGER-DECODE-FLOAT die
  // drei Zahlen m,e,s (Mantisse, Exponent, Vorzeichen) liefert,
  // = s * 2^e * m.
  // n:=m. Falls s<0, setze n:=-m.
  // Falls e>=0, ist (ash n e) das Ergebnis,
  // sonst ist die rationale Zahl (/ n (ash 1 (- e))) das Ergebnis.
	var cl_idecoded_float x_decoded = integer_decode_float(x);
	var cl_I& m = x_decoded.mantissa;
	var cl_I& e = x_decoded.exponent;
	var cl_I& s = x_decoded.sign;
	var cl_I n = (!minusp(s) ? m : -m);
	if (!minusp(e))
		return ash(n,e);
	else {
#if 0
		return I_posI_div_RA(n, ash(1,-e));
#else // spart ggT
		// n /= 0, -e > 0. Kürze mit ggT(n,2^(-e)) = 2^min(ord2(n),-e).
		// 0 < -e <= LF_exp_mid-LF_exp_low + intDsize*len < 2^32,
		var cl_I minus_e = -e;
		var uintL _e = cl_I_to_UL(minus_e); // daher kein Überlauf
		var uintC k = ord2(n);
		if (k >= _e)
			// Kürze mit 2^(-e).
			return ash(n,e);
		else
			// Kürze mit 2^k, 0 <= k < -e.
			return I_I_to_RT(ash(n,-(sintC)k),
					 ash(1,minus_e-(cl_I)(unsigned long)k));
#endif
	}
}

}  // namespace cln
