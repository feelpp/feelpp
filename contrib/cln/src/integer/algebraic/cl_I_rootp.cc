// rootp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

// Methode:
// Falls x=0 oder x=1: x = x^n -> JA, x als Ergebnis.
// Hier also x>1. Suche ein Integer y > 1 mit x=y^n.
// Falls n >= integer_length(x): NEIN. (Da y>=2, müßte x>=2^n gelten.)
// Hier also n>0 klein...

bool rootp (const cl_I& x, uintL n, cl_I* w)
{
	if (eq(x,0) || eq(x,1)) // x=0 oder x=1 ?
	  { *w = x; return true; } // ja -> x als Ergebnis
	if (n >= integer_length(x))
	  { return false; }
	return cl_rootp_aux(x,n,w);
}

}  // namespace cln
