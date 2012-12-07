// div2adic().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digit/cl_2D.h"


// Implementation.

namespace cln {

uintD div2adic (uintD a, uintD b)
{
// Methode:
// Konstruiere c Bit für Bit.
// c := 0, d := a.
// Für j=0,...,intDsize:
//   [Hier b*c == a mod 2^j und d = (a-b*c)/2^j.] j=intDsize -> fertig.
//   Falls d ungerade, setze c:=c+2^j und d:=(d-b)/2, sonst d:=d/2.
// Ergebnis c.
      ASSERT(!((b % 2) ==0))
#if 1
     {var uintD c = 0;
      var uintD bit_j = 1; // 2^j
      loop // Verwende a als Variable d
        { if (a & bit(0)) { c = c+bit_j; a = a-b; }
          a = a>>1;
          bit_j = bit_j << 1;
          if (bit_j == 0) break; // j=intDsize -> fertig
        }
      return c;
     }
#else
     {var uintD bit_j = 1; // 2^j
      var uintD b_j = b-1; // (b-1)*2^j
      loop // Verwende a als Variable d*2^j+c
        { if (a & bit_j) { a = a - b_j; }
          b_j = b_j << 1; bit_j = bit_j << 1;
          if (bit_j == 0) break; // j=intDsize -> fertig
        }
      return a;
     }
#endif
}

}  // namespace cln
