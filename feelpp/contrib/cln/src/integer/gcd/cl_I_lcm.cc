// lcm().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I lcm (const cl_I& a, const cl_I& b)
{
// Methode:
// a=0 oder b=0 -> Ergebnis 0.
// a:=(abs a), b:=(abs b).
// g:=ggT(a,b)>0.
// Falls g=1, Ergebnis a*b, sonst Ergebnis (a/g)*b.
      if (eq(a,0)) { return 0; }
      if (eq(b,0)) { return 0; }
      // Beträge nehmen:
     {var cl_I abs_a = abs(a);
      var cl_I abs_b = abs(b);
      var cl_I& a = abs_a;
      var cl_I& b = abs_b;
      var cl_I g = gcd(a,b);
      if (!eq(g,1)) { a = exquopos(a,g); } // a durch g (beide >0) dividieren
      return a*b;
     }
}

}  // namespace cln
