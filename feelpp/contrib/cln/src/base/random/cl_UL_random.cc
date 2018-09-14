// Word level random number generator.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/random.h"


// Implementation.

#include "base/cl_low.h"

namespace cln {

// Zufallszahlengenerator nach [Knuth: The Art of Computer Programming, Vol. II,
// Seminumerical Algorithms, 3.3.4., Table 1, Line 30], nach C. Haynes:
// X eine 64-Bit-Zahl. Iteration X := (a*X+c) mod m
// mit m=2^64, a=6364136223846793005, c=1.

uint32 random32 (random_state& randomstate)
{
#ifdef HAVE_FAST_LONGLONG
	// Multiplikator a=6364136223846793005 = 0x5851F42D4C957F2D :
	var uint64 seed = highlow64(randomstate.seed.hi,randomstate.seed.lo);
	var const uint64 a = 0x5851F42D4C957F2DULL;
	var uint64 newseed;
	// multiplizieren, brauche nur letzte 64 Bit:
	mulu64(seed,a, , newseed =);
	// c addieren:
	newseed += 1;
	// seed neu füllen:
	randomstate.seed.hi = high32(newseed);
	randomstate.seed.lo = low32(newseed);
	// mittlere 32 Bits als Ergebnis:
	return (uint32)(newseed >> 16);

#else
	// Multiplikator a=6364136223846793005 = 0x5851F42D4C957F2D :
	var uint32 seed_hi = randomstate.seed.hi;
	var uint32 seed_lo = randomstate.seed.lo;
	var const uint32 a_hi = 0x5851F42D;
	var const uint32 a_lo = 0x4C957F2D;
	var uint32 newseed_hi;
	var uint32 newseed_lo;
	// multiplizieren, brauche nur letzte 64 Bit:
	mulu32(seed_lo,a_lo, newseed_hi =, newseed_lo =);
	mulu32(seed_lo,a_hi, , newseed_hi +=);
	mulu32(seed_hi,a_lo, , newseed_hi +=);
	// c addieren:
	newseed_lo += 1; if (newseed_lo==0) { newseed_hi += 1; } // um 1 erhöhen
	// seed neu füllen:
	randomstate.seed.hi = newseed_hi;
	randomstate.seed.lo = newseed_lo;
	// mittlere 32 Bits als Ergebnis:
	return highlow32(low16(newseed_hi),high16(newseed_lo));
#endif
}

}  // namespace cln
