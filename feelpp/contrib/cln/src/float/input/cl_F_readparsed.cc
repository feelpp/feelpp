// read_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float_io.h"


// Implementation.

#include "cln/integer.h"
#include "integer/cl_I.h"
#include "cln/rational.h"
#include "rational/cl_RA.h"
#include "cln/float.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"
#include "float/cl_F.h"
#include "float/sfloat/cl_SF.h"
#include "float/ffloat/cl_FF.h"
#include "float/dfloat/cl_DF.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F read_float (unsigned int base, float_format_t prec, cl_signean sign, const char * string, uintC index1, uintC index4, uintC index2, uintC index3)
{
	var cl_I exponent;
	{
		var uintC exp_len = index2-index4; // Anzahl Stellen des Exponenten
		if (exp_len > 0) {
			var const char * ptr = &string[index4]; // zeigt auf den Exponentmarker
			ptr++; exp_len--; // Exponentmarker überlesen
			var cl_signean exp_sign = 0; // Exponenten-Vorzeichen
			switch (*ptr) {
				case '-': exp_sign = ~exp_sign; // fallthrough
				case '+': ptr++; exp_len--; // Exponenten-Vorzeichen überlesen
				default: ;
			}
			exponent = digits_to_I(ptr,exp_len,(uintD)base); // Exponent in Integer umwandeln
			if (!(exp_sign==0))
				exponent = -exponent; // incl. Vorzeichen
		} else {
			// kein Exponent da
			exponent = 0;
		}
	}
	// exponent = Exponent als Integer.
	var cl_RA base_power = expt(base,exponent-(index4-index3)); // zu multiplizierende Zehnerpotenz
	var cl_I mantisse = // Mantisse als Integer
	  digits_to_I(&string[index1],index4-index1,(uintD)base);
	// Mantisse (Integer) und Zehnerpotenz (rational >0) unelegant zusammenmultiplizieren:
	var cl_RA ratvalue;
	if (integerp(base_power)) {
		DeclareType(cl_I,base_power);
		ratvalue = mantisse * base_power;
	} else {
		// falls mantisse/=0, in exponent=1/10^i den Zähler durch mantisse
		// ersetzen (liefert ungekürzten Bruch, Vorsicht!)
		DeclareType(cl_RT,base_power);
		if (zerop(mantisse))
			ratvalue = 0;
		else {
			ASSERT(TheRatio(base_power)->refcount == 1);
			TheRatio(base_power)->numerator = mantisse;
			ratvalue = base_power;
		}
	}
	// ratvalue = Mantisse * Zehnerpotenz, als ungekürzte rationale Zahl!
	floatformatcase(prec
	,	// in Short-Float umwandeln
		{
			var cl_SF x = cl_RA_to_SF(ratvalue);
			if (sign==0) { return x; } else { return -x; }
		}
	,	// in Single-Float umwandeln
		{
			var cl_FF x = cl_RA_to_FF(ratvalue);
			if (sign==0) { return x; } else { return -x; }
		}
	,	// in Double-Float umwandeln
		{
			var cl_DF x = cl_RA_to_DF(ratvalue);
			if (sign==0) { return x; } else { return -x; }
		}
	,	// in Long-Float umwandeln
		{
			var cl_LF x = cl_RA_to_LF(ratvalue,len);
			if (sign==0) { return x; } else { return -x; }
		}
	);
}

}  // namespace cln
