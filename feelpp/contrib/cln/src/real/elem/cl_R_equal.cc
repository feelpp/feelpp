// equal().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "rational/cl_RA.h"
#include "cln/integer.h"

namespace cln {

inline bool equal (const cl_F& x, const cl_F& y)
{
	return compare(x,y) == 0;
}

bool equal (const cl_R& x, const cl_R& y)
{
// Methode:
// Beide rational oder beide Floats -> klar.
// Eine rational, eine Float ->
//   Die rationale Zahl muß einen Zweierpotenz-Nenner haben, sonst verschieden.
//   Die rationale Zahl zum Float machen, vergleichen.
//   Verschieden -> Das war's.
//   Gleich -> Das Float mit RATIONAL rational machen, nochmals vergleichen.
	realcase2(x
	,	realcase2(y
		,	// beides rationale Zahlen
			return equal(x,y);
		,	// x rational, y Float -> x in Float umwandeln
			if (!power2p(denominator(x)))
				return false;
			if (!equal(cl_float(x,y),y))
				return false;
			return equal(x,rational(y));
		);
	,	realcase2(y
		,	// x Float, y rational -> y in Float umwandeln
			if (!power2p(denominator(y)))
				return false;
			if (!equal(x,cl_float(y,x)))
				return false;
			return equal(rational(x),y);
		,	// beides Floats
			return equal(x,y);
		);
	);
}

}  // namespace cln
