// expt_pos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digit/cl_2D.h"


// Implementation.

namespace cln {

uintD expt_pos (uintD a, uintL b)
{
// Methode:
//   a:=x, b:=y, c:=1. [a^b*c bleibt invariant, = x^y.]
//   Solange b>1,
//     falls b ungerade, setze c:=a*c,
//     setze b:=floor(b/2),
//     setze a:=a*a.
//   Wenn b=1, setze c:=a*c.
//   Liefere c.
// Oder optimiert:
//   a:=x, b:=y.
//   Solange b gerade, setze a:=a*a, b:=b/2. [a^b bleibt invariant, = x^y.]
//   c:=a.
//   Solange b:=floor(b/2) >0 ist,
//     setze a:=a*a, und falls b ungerade, setze c:=a*c.
//   Liefere c.
      while ((b & bit(0)) == 0) { a = mul2adic(a,a); b = b>>1; }
      var uintD c = a;
      until ((b = b>>1) == 0)
        { a = mul2adic(a,a);
          if (b & bit(0)) { c = mul2adic(a,c); }
        }
      return c;
}

}  // namespace cln
