// float_sign().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

namespace cln {

const cl_F float_sign (const cl_F& x, const cl_F& y)
{
  // Methode:
  // Falls x<0 xor y<0, Ergebnis (- y), sonst Ergebnis y.
	if (minusp(x) != minusp(y))
		return -y;
	else
		return y;
}

}  // namespace cln
