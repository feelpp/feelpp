// cornacchia4().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/numtheory.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

// [Cohen], section 1.5.2, algorithm 1.5.3.
// For proofs refer to [F. Morain, J.-L. Nicolas: On Cornacchia's algorithm
// for solving the diophantine equation u^2+v*d^2=m].

const cornacchia_t cornacchia4 (const cl_I& d, const cl_I& p)
{
	// Method:
	// Goal: Solve x^2+d*y^2 = 4*p.
	// If p=2: {0,1,4,...} + d*{0,1,4,...} = 8.
	//   d=1: (x,y) = (2,2).
	//   d=2: (x,y) = (0,2).
	//   d=4: (x,y) = (2,1).
	//   d=7: (x,y) = (1,1).
	//   Else no solution.
	// If p>2:
	//   If d == 0 mod 4:
	//     x must be even. Solve u^2 + (d/4)*y^2 = p, return (2*u,y).
	//   If d == 2 mod 4:
	//     x must be even, then y must be even. Solve u^2 + d*v^2 = p,
	//     return (2*u,2*v).
	//   If d == 1,5,7 mod 8:
	//     x^2+d*y^2 == 4 mod 8, x^2 and y^2 can only be == 0,1,4 mod 8.
	//     x and y must be even. Solve u^2 + d*v^2 = p, return (2*u,2*v).
	//   If d == 3 mod 8:
	//     Compute x with x^2+d == 0 mod 4*p.
	//     Euclidean algorithm on a = 2*p, b = x, stop when a remainder
	//     <= 2*sqrt(p) has been reached.
	var cl_I p4 = p<<2;
	if (d >= p4) {
		if (d == p4)
			// (x,y) = (0,1)
			return cornacchia_t(1, 0,1);
		else
			// d > 4*p -> no solution
			return cornacchia_t(0);
	}
	// Now 0 < d < 4*p.
	if (p == 2) {
		if (d==1) return cornacchia_t(1, 2,2);
		if (d==2) return cornacchia_t(1, 0,2);
		if (d==4) return cornacchia_t(1, 2,1);
		if (d==7) return cornacchia_t(1, 1,1);
		return cornacchia_t(0);
	}
	switch (FN_to_V(logand(d,7))) {
		case 0: case 4: {
			// d == 0 mod 4
			var cornacchia_t s = cornacchia1(d>>2,p);
			if (!s.condition)
				if (s.solutions != 0)
					s.solution_x = s.solution_x<<1;
			return s;
		}
		case 1: case 2: case 5: case 6: case 7: {
			var cornacchia_t s = cornacchia1(d,p);
			if (!s.condition)
				if (s.solutions != 0) {
					s.solution_x = s.solution_x<<1;
					s.solution_y = s.solution_y<<1;
				}
			return s;
		}
		case 3:
			break;
	}
	switch (jacobi(-d,p)) {
		case -1: // no solution
			return cornacchia_t(0);
		case 0: // gcd(d,p) > 1
			return new cl_composite_condition(p,gcd(d,p));
		case 1:
			break;
	}
	// Compute x with x^2+d == 0 mod p.
	var cl_modint_ring R = find_modint_ring(p);
	var sqrt_mod_p_t init = sqrt_mod_p(R,R->canonhom(-d));
	if (init.condition)
		return init.condition;
	if (init.solutions != 2)
		throw runtime_exception();
	// Compute x with x^2+d == 0 mod 4*p.
	var cl_I x0 = R->retract(init.solution[0]);
	if (evenp(x0)) { x0 = p-x0; } // Enforce x0^2+d == 0 mod 4.
	// Euclidean algorithm.
	var cl_I a = p<<1;
	var cl_I b = x0;
	var cl_I limit = isqrt(p4);
	while (b > limit) {
		var cl_I r = mod(a,b);
		a = b; b = r;
	}
	// b is the first euclidean remainder <= 2*sqrt(p).
	var cl_I& x = b;
	var cl_I_div_t div = floor2(p4-square(b),d);
	if (!zerop(div.remainder))
		return cornacchia_t(0);
	var cl_I& c = div.quotient;
	var cl_I y;
	if (!sqrtp(c,&y))
		return cornacchia_t(0);
	return cornacchia_t(1, x,y);
}

}  // namespace cln
