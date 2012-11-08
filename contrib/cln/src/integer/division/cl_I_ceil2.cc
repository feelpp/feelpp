// ceiling2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I_div_t ceiling2 (const cl_I& x, const cl_I& y)
{
// Methode:
// (ceiling x y) :==
// (DIVIDE (abs x) (abs y)) -> q,r
// Falls x,y selbes Vorzeichen haben und r<>0,
//   setze q:=q+1 und r:=r-abs(y).
// Falls x<0, setze r:=-r.
// Falls x,y verschiedene Vorzeichen haben, setze q:=-q.
// Liefere q,r.
  var cl_I abs_y = abs(y);
  var cl_I_div_t q_r = cl_divide(abs(x),abs_y);
  var cl_I& q = q_r.quotient;
  var cl_I& r = q_r.remainder;
  if ((minusp(x) == minusp(y)) && !zerop(r))
    { q = q + 1; r = r - abs_y; }
  if (minusp(x))
    { r = -r; }
  if (minusp(x) != minusp(y))
    { q = -q; }
  return q_r;
}

}  // namespace cln
