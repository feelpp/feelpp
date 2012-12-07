// Digit level 2-adic arithmetic

#ifndef _CL_2D_H
#define _CL_2D_H

#include "cln/types.h"
#include "base/digit/cl_D.h"

namespace cln {

// Multipliziert zwei Zahlen mod 2^intDsize.
// mul2adic(a,b)
// > uintD a,b: Zahlen mod 2^intDsize
// < ergebnis: Zahl c mod 2^intDsize mit c == a*b mod 2^intDsize
  extern uintD mul2adic (uintD a, uintD b);
#if HAVE_DD
  inline uintD mul2adic (uintD a, uintD b)
  {
	return lowD(muluD(a,b));
  }
#else
  inline uintD mul2adic (uintD a, uintD b)
  {
	muluD(a,b, ,return);
  }
#endif

// Potenziert eine Zahl mod 2^intDsize.
// expt_pos(x,y)
// > uintD x: Zahl mod 2^intDsize
// > uintL y: Exponent >0
// < uintD ergebnis: x^y mod 2^intDsize
  extern uintD expt_pos (uintD x, uintL y);

// Dividiert zwei Zahlen mod 2^intDsize.
// div2adic(a,b)
// > uintD a: Zahl mod 2^intDsize
// > uintD b: ungerade Zahl mod 2^intDsize
// < ergebnis: Zahl c mod 2^intDsize mit b*c == a mod 2^intDsize
  extern uintD div2adic (uintD a, uintD b);

}  // namespace cln

#endif /* _CL_2D_H */
