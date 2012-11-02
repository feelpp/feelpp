// cl_LF_len_incsqrtx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

namespace cln {

uintC cl_LF_len_incsqrtx (uintC n)
{
// Methode bei intDsize=16:
// Allgemein: intDsize*n + sqrt(intDsize*n) + 2 + 31 < intDsize*(n+inc)
// <==>       sqrt(intDsize*n) + 33 < intDsize*inc
// <==>       sqrt(intDsize*n) < intDsize*inc - 33
// <==>       intDsize*n < intDsize^2*inc^2 - 66*intDsize*inc + 1089
// <==>       n <= intDsize*inc^2 - 66*inc + floor(1089/intDsize)
	return
	  #define NMAX(k)  (uintC)((intDsize*(k)-66)*(k)+floor(1089,intDsize))
	  #define FITS(n,k)  ((intDsize*(k) > 33) && ((n) <= NMAX(k)))
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
