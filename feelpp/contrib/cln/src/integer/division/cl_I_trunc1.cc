// truncate1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I truncate1 (const cl_I& x, const cl_I& y)
{
// Methode:
// (truncate x y) :==
// (DIVIDE (abs x) (abs y)) -> q,r
// Falls x<0, setze r:=-r.
// Falls x,y verschiedene Vorzeichen haben, setze q:=-q.
// Liefere nur q.
  var cl_I_div_t q_r = cl_divide(abs(x),abs(y));
  var cl_I& q = q_r.quotient;
  if (minusp(x) != minusp(y))
    { q = -q; }
  return q;
}

}  // namespace cln
