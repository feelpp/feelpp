// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "cln/float.h"

namespace cln {

cl_signean compare (const cl_R& x, const cl_R& y)
{
  // Methode:
  // Beide rational oder beide Floats -> klar.
  // Eine rational, eine Float ->
  //   Die rationale Zahl zum Float machen, vergleichen.
  //   Verschieden -> Das war's.
  //   Gleich -> Das Float mit RATIONAL rational machen, nochmals vergleichen.
	realcase2(x
	,	realcase2(y
		,	// beides rationale Zahlen
			return compare(x,y);
		,	// x rational, y Float -> x in Float umwandeln
			var cl_signean result = compare(cl_float(x,y),y);
			if (result != signean_null)
				return result;
			return compare(x,rational(y));
		);
	,	realcase2(y
		,	// x Float, y rational -> y in Float umwandeln
			var cl_signean result = compare(x,cl_float(y,x));
			if (result != signean_null)
				return result;
			return compare(rational(x),y);
		,	// beides Floats
			return compare(x,y);
		);
	);
}

}  // namespace cln
