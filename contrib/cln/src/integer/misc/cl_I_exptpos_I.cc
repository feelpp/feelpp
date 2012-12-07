// expt_pos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I expt_pos (const cl_I& x, const cl_I& y)
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
  #if 0 // unoptimiert
    var cl_I a = x;
    var cl_I b = y;
    var cl_I c = 1;
    until (eq(b,1))
      { if (oddp(b))
          { c = a * c; }
        b = b >> 1; // b := (floor b 2)
        a = square(a);
      }
    return a * c;
  #else // optimiert
    var cl_I a = x;
    var cl_I b = y;
    while (!oddp(b)) { a = square(a); b = b >> 1; }
    var cl_I c = a;
    until (eq(b,1))
      { b = b >> 1;
        a = square(a);
        if (oddp(b)) { c = a * c; }
      }
    return c;
  #endif
}
// Bit complexity (x of length N): O(M(N*y)).

}  // namespace cln
