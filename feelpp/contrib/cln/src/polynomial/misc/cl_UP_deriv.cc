// deriv().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly.h"


// Implementation.

#include "cln/integer.h"

namespace cln {

const cl_UP deriv (const cl_UP& x)
{
	// Method:
	// Write x = a0 T^0 + ... + an T^n.
	// Then deriv(x) = 1*a1 T^0 + ... + n*an T^(n-1)  (= 0 if n <= 0).
	var cl_univpoly_ring UPR = x.ring();
	var sintL n = degree(x);
	if (n <= 0)
		return UPR->zero();
	else {
		var cl_UP y = UPR->create(n-1);
		for ( ; n > 0; n--)
			y.set_coeff(n-1, n * coeff(x,n));
		y.finalize();
		return y;
	}
}

}  // namespace cln
