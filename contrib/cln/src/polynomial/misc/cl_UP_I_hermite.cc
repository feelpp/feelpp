// hermite().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly_integer.h"


// Implementation.

#include "cln/integer.h"

namespace cln {

const cl_UP_I hermite (sintL n)
{
// The Hermite polynomials H_n(x) are defined as
//
//                             ( d  ) n
//    H_n(x) = (-1)^n exp(x^2) (----)   exp(- x^2)
//                             ( dx )
//
// They satisfy the recurrence relation
//
//    H_0(x) = 1
//    H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x) for n >= 0.
//
// Theorem:
//    H_n(x) satisfies the differential equation
//    H_n''(x) - 2x*H_n'(x) + 2n*H_n(x) = 0.
//
// Proof: See elsewhere.
//
// Corollary:
//    The coefficients c_{n,k} of H_n(x) = sum(k=0..n, c_{n,k} x^k)
//    satisfy:
//       c_{n,n} = 2^n,
//       c_{n,n-1} = 0,
//       c_{n,k} = (k+1)(k+2)/2(k-n)*c_{n,k+2}
//
// It follows that for n>=0
//
//    H_n(x) = sum(j=0..floor(n/2), (-1)^j n!/j!(n-2j)! 2^(n-2j) x^(n-2j))
//
	var cl_univpoly_integer_ring R = find_univpoly_ring(cl_I_ring);
	var cl_UP_I h = R->create(n);
	var sintL k = n;
	var cl_I c_k = ash(1,n);
	for (;;) {
		h.set_coeff(k,c_k);
		k = k-2;
		if (k < 0)
			break;
		c_k = exquo((cl_I)(k+1) * (cl_I)(k+2) * c_k,
		            2*(cl_I)(k-n));
	}
	h.finalize();
	return h;
}

}  // namespace cln
