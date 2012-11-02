// round1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I round1 (const cl_I& x, const cl_I& y)
{
// Methode:
// (round x y) :==
// (DIVIDE (abs x) (abs y)) -> q,r
// Setze s:=abs(y)-r.
// Falls (r>s) oder (r=s und q ungerade),
//   (d.h. falls r>abs(y)/2 oder r=abs(y)/2 und q ungerade),
//   setze q:=q+1 und r:=-s (d.h. r:=r-abs(y)).
// {Nun ist abs(r) <= abs(y)/2, bei abs(r)=abs(y)/2 ist q gerade.}
// Falls x<0, setze r:=-r.
// Falls x,y verschiedene Vorzeichen haben, setze q:=-q.
// Liefere nur q.
  var cl_I abs_y = abs(y);
  var cl_I_div_t q_r = cl_divide(abs(x),abs_y);
  var cl_I& q = q_r.quotient;
  var cl_I& r = q_r.remainder;
  var cl_I s = abs_y - r;
  if ((r > s) || ((r == s) && oddp(q)))
    { q = q + 1; }
  if (minusp(x) != minusp(y))
    { q = -q; }
  return q;
}

}  // namespace cln
