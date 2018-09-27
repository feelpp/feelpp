// exquopos().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I exquopos (const cl_I& x, const cl_I& y)
{
// Methode:
// (exquopos x y) :==
// (DIVIDE x y) -> q,r
// Falls r<>0, Error.
// Liefere q.
  var cl_I_div_t q_r = cl_divide(x,y);
  if (!zerop(q_r.remainder)) { throw exquo_exception(x,y); }
  return q_r.quotient;
}

}  // namespace cln
