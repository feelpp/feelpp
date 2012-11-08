// cl_LF_len_incsqrt().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

namespace cln {

uintC cl_LF_len_incsqrt (uintC n)
{
// Methode bei intDsize=16:
// n -> n+1 für n<=12 wegen 16n+sqrt(16n)+2 < 16(n+1)
// n -> n+2 für n<=56 wegen 16n+sqrt(16n)+2 < 16(n+2)
// n -> n+4 für n<=240
// n -> n+8 für n<=992
// n -> n+16 für n<=4032
// n -> n+32 für n<=16256
// n -> n+65 für n<=65535
// Allgemein: intDsize*n + sqrt(intDsize*n) + 2 < intDsize*(n+inc)
// <==>       sqrt(intDsize*n) + 2 < intDsize*inc
// <==>       sqrt(intDsize*n) < intDsize*inc - 2
// <==>       intDsize*n < intDsize^2*inc^2 - 4*intDsize*inc + 4
// <==>       n <= intDsize*inc^2 - 4*inc
	return
	  #define NMAX(k)  (uintC)((intDsize*(k)-4)*(k))
	  #define FITS(n,k)  ((n) <= NMAX(k))
	  #define n_max  (uintC)(bitm(intCsize)-1)
	  #define TRYINC(inc)  FITS(n_max,inc) || FITS(n,inc) ? RETINC(inc) :
	  #define RETINC(inc)  \
	    /* at this point we know  n <= NMAX(inc) */			\
	    /* Check whether n + (inc) overflows. */			\
	  ((NMAX(inc) <= n_max-(inc)) || (n <= n_max-(inc)) ? n+(inc) : n_max)
	  #define TEST(i)  TRYINC(1UL<<i)
	  TEST(0)
	  TEST(1)
	  TEST(2)
	  TEST(3)
	  TEST(4)
	  TEST(5)
	  TEST(6)
	  TEST(7)
	  TEST(8)
	  TEST(9)
	  TEST(10)
	  TEST(11)
	  TEST(12)
	  TEST(13)
	  // No TEST(14), because NMAX(1UL<<14) is already out of range.
	  n_max;
}

}  // namespace cln
