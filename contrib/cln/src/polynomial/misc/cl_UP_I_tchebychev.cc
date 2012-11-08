// tschebychev().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly_integer.h"


// Implementation.

#include "cln/integer.h"

namespace cln {

const cl_UP_I tschebychev (sintL n)
{
// The Tschebychev polynomials (of the 1st kind) T_n(x) are defined
// through the recurrence relation
//
//    T_0(x) = 1
//    T_1(x) = x
//    T_{n+2}(x) = 2x T_{n+1}(x) - T_n(x) for n >= 0.
//
// Theorem:
//    T_n(x) satisfies the differential equation
//    (x^2-1)*T_n''(x) + x*T_n'(x) - n^2*T_n(x) = 0.
//
// Proof: See elsewhere.
//
// Corollary:
//    The coefficients c_{n,k} of T_n(x) = sum(k=0..n, c_{n,k} x^k)
//    satisfy:
//       c_{n,n} = 2^(n-1) for n>=1, 1 for n=0,
//       c_{n,n-1} = 0,
//       c_{n,k} = (k+1)(k+2)/(k^2-n^2)*c_{n,k+2}
//
// It follows that for n>0
//
// T_n(x) = sum(j=0..floor(n/2), (-1)^j (n-j-1)!n/j!(n-2j)! 2^(n-2j-1) x^(n-2j))
//
	var cl_univpoly_integer_ring R = find_univpoly_ring(cl_I_ring);
	if (n == 0)
		return R->one();
	var cl_UP_I t = R->create(n);
	var sintL k = n;
	var cl_I c_k = ash(1,n-1);
	for (;;) {
		t.set_coeff(k,c_k);
		k = k-2;
		if (k < 0)
			break;
		c_k = exquo((cl_I)(k+1) * (cl_I)(k+2) * c_k,
		            (cl_I)(k-n) * (cl_I)(k+n));
	}
	t.finalize();
	return t;
}

}  // namespace cln
