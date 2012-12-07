// legendre().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/univpoly_rational.h"


// Implementation.

#include "cln/integer.h"
#include "cln/rational.h"

namespace cln {

const cl_UP_RA legendre (sintL n)
{
// The Legendre polynomials P_n(x) are defined as
//
//                1   ( d  ) n          n
//    P_n(x) = ------ (----)   (x^2 - 1)
//             2^n n! ( dx )
//
// They satisfy the recurrence relation
//
//    P_0(x) = 1
//    P_{n+1}(x) = 1/(n+1) * ((2n+1)x P_n(x) - n P_{n-1}(x)) for n >= 0.
//
// Theorem:
//    P_n(x) satisfies the differential equation
//    (1-x^2)*P_n''(x) - 2x*P_n'(x) + (n^2+n)*P_n(x) = 0.
//
// Proof: See elsewhere.
//
// Corollary:
//    The coefficients c_{n,k} of P_n(x) = sum(k=0..n, c_{n,k} x^k)
//    satisfy:
//       c_{n,n} = (2n)!/(n!^2 2^n),
//       c_{n,n-1} = 0,
//       c_{n,k} = (k+1)(k+2)/(k-n)(k+n+1)*c_{n,k+2}
//
// It follows that for n>=0
//
// P_n(x) = sum(j=0..floor(n/2), (-1)^j (2n-2j)!/j!(n-2j)!(n-j)! 2^-n x^(n-2j))
//
	var cl_univpoly_rational_ring R = find_univpoly_ring(cl_RA_ring);
	var cl_UP_RA p = R->create(n);
	var cl_I denom = ash(1,n);
	var sintL k = n;
	var cl_I c_k = binomial(2*n,n);
	for (;;) {
		p.set_coeff(k,c_k/denom);
		k = k-2;
		if (k < 0)
			break;
		c_k = exquo((cl_I)(k+1) * (cl_I)(k+2) * c_k,
		            (cl_I)(k-n) * (cl_I)(k+n+1));
	}
	p.finalize();
	return p;
}

}  // namespace cln
