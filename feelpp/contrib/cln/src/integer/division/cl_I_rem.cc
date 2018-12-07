// rem().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I rem (const cl_I& x, const cl_I& y)
{
// Methode:
// (rem x y) :==
// (DIVIDE (abs x) (abs y)) -> q,r
// Falls x<0, setze r:=-r.
// Liefere r.
  var cl_I r = cl_divide(abs(x),abs(y)).remainder;
  if (minusp(x)) { return -r; } else { return r; }
}

}  // namespace cln
