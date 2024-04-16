// float_format().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

namespace cln {

float_format_t float_format (uintE n)
{
// Methode:
// Mindestens 1+n Dezimalstellen (inklusive Vorkommastelle)
// bedeutet mindestens ceiling((1+n)*ln(10)/ln(2)) Binärstellen.
// ln(10)/ln(2) = 3.321928095 = (binär) 11.010100100110100111100001001011110011010...
//                       = (binär) 100 - 0.101011011001011000011110110100001100101...
// Durch diese Berechnungsmethode wird das Ergebnis sicher >= (1+n)*ln(10)/ln(2)
// sein, evtl. um ein paar Bit zu groß aber nicht zu klein.
	n = 1+n;
	return (float_format_t)
	       ((n << 2)
		- (n >> 1) - (n >> 3) - (n >> 5) - (n >> 6)
		- (n >> 8) - (n >> 9) - (n >> 12) - (n >> 14)
		- (n >> 15) - (n >> 20) - (n >> 21) - (n >> 22)
		- (n >> 23) - (n >> 25) - (n >> 26) - (n >> 28)
#if (intEsize>32)
		- (n >> 33) - (n >> 34) - (n >> 37) - (n >> 39)
#endif
		);
}

}  // namespace cln
