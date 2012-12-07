// mod().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I mod (const cl_I& x, const cl_I& y)
{
// Methode:
// (mod x y) :==
// (DIVIDE (abs x) (abs y)) -> q,r
// Falls x,y verschiedene Vorzeichen haben und r<>0, setze r:=r-abs(y).
// Falls x<0, setze r:=-r.
// Liefere r.
  var cl_I abs_y = abs(y);
  var cl_I r = cl_divide(abs(x),abs_y).remainder;
  if (minusp(x) != minusp(y))
    { if (zerop(r)) { return 0; }
      r = r - abs_y;
    }
  if (minusp(x)) { return -r; } else { return r; }
}

}  // namespace cln
