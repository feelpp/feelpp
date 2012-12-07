// laguerre().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly_integer.h"


// Implementation.

#include "cln/integer.h"

namespace cln {

const cl_UP_I laguerre (sintL n)
{
// The Laguerre polynomials L_n(x) are defined as
//
//                    ( d  ) n
//    L_n(x) = exp(x) (----)   (x^n exp(-x))
//                    ( dx )
//
// They satisfy the recurrence relation
//
//    L_0(x) = 1
//    L_{n+1}(x) = (2n+1-x) L_n(x) - n^2 L_{n-1}(x) for n >= 0.
//
// Theorem:
//    L_n(x) satisfies the differential equation
//    x*L_n''(x) + (1-x)*L_n'(x) + n*L_n(x) = 0.
//
// Proof: See elsewhere.
//
// Corollary:
//    The coefficients c_{n,k} of L_n(x) = sum(k=0..n, c_{n,k} x^k)
//    satisfy:
//       c_{n,n} = (-1)^n,
//       c_{n,k} = (k+1)^2/(k-n)*c_{n,k+1}
//
// It follows that for n>=0
//
//       L_n(x) = sum(j=0..n, (-1)^(n-j) n!^2/j!(n-j)!^2 x^(n-j))
//
	var cl_univpoly_integer_ring R = find_univpoly_ring(cl_I_ring);
	var cl_UP_I l = R->create(n);
	var sintL k = n;
	var cl_I c_k = (evenp(n) ? 1 : -1);
	for (;;) {
		l.set_coeff(k,c_k);
		k = k-1;
		if (k < 0)
			break;
		c_k = exquo((cl_I)(k+1) * (cl_I)(k+1) * c_k,
		            (cl_I)(k-n));
	}
	l.finalize();
	return l;
}

}  // namespace cln
