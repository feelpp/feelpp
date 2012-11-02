// exquo().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I exquo (const cl_I& x, const cl_I& y)
{
// Methode:
// (exquo x y) :==
// (DIVIDE (abs x) (abs y)) -> q,r
// Falls r<>0, Error.
// Falls x,y verschiedene Vorzeichen haben, liefere -q, sonst q.
  var cl_I_div_t q_r = cl_divide(abs(x),abs(y));
  if (!zerop(q_r.remainder)) { throw exquo_exception(x,y); }
  if (minusp(x) == minusp(y))
    { return q_r.quotient; }
    else
    { return -q_r.quotient; }
}

}  // namespace cln
