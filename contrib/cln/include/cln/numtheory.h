// Number theoretic operations.

#ifndef _CL_NUMTHEORY_H
#define _CL_NUMTHEORY_H

#include "cln/number.h"
#include "cln/integer.h"
#include "cln/modinteger.h"
#include "cln/condition.h"

namespace cln {

// jacobi(a,b) returns the Jacobi symbol
//                       (  a  )
//                       ( --- )
//                       (  b  )
// a, b must be integers, b > 0, b odd. The result is 0 iff gcd(a,b) > 1.
  extern int jacobi (sintV a, sintV b);
  extern int jacobi (const cl_I& a, const cl_I& b);

// isprobprime(n), n integer > 0,
// returns true when n is probably prime.
// This is pretty quick, but no caching is done.
  extern bool isprobprime (const cl_I& n);

// nextprobprime(x) returns the smallest probable prime >= x.
  extern const cl_I nextprobprime (const cl_R& x);

#if 0
// primitive_root(R) of R = Z/pZ, with p a probable prime,
// returns
//   either a generator of (Z/pZ)^*, assuming p is prime, or
//   a proof that p is not prime, maybe even a non-trivial factor of p.
struct primitive_root_t {
	cl_composite_condition* condition;
	cl_MI gen;
	// Constructors.
	primitive_root_t (cl_composite_condition* c) : condition (c) {}
	primitive_root_t (const cl_MI& g) : condition (NULL), gen (g) {}
};
extern const primitive_root_t primitive_root (const cl_modint_ring& R);
#endif

// sqrt_mod_p(R,x) where x is an element of R = Z/pZ, with p a probable prime,
// returns
//   either the square roots of x in R, assuming p is prime, or
//   a proof that p is not prime, maybe even a non-trivial factor of p.
struct sqrt_mod_p_t {
	cl_composite_condition* condition;
	// If no condition:
	int solutions; // 0,1,2
	cl_I factor; // zero or non-trivial factor of p
	cl_MI solution[2]; // max. 2 solutions
	// Constructors.
	sqrt_mod_p_t () {}
	sqrt_mod_p_t (cl_composite_condition* c) : condition (c) {}
	sqrt_mod_p_t (int s) : condition (NULL), solutions (s) {}
	sqrt_mod_p_t (int s, const cl_MI& x0) : condition (NULL), solutions (s)
		{ solution[0] = x0; }
	sqrt_mod_p_t (int s, const cl_MI& x0, const cl_MI& x1) : condition (NULL), solutions (s)
		{ solution[0] = x0; solution[1] = x1; }
};
extern const sqrt_mod_p_t sqrt_mod_p (const cl_modint_ring& R, const cl_MI& x);

// cornacchia1(d,p) solves x^2 + d*y^2 = p.
// cornacchia4(d,p) solves x^2 + d*y^2 = 4*p.
// d is an integer > 0, p is a probable prime.
// It returns
//   either a nonnegative solution (x,y), if it exists, assuming p is prime, or
//   a proof that p is not prime, maybe even a non-trivial factor of p.
struct cornacchia_t {
	cl_composite_condition* condition;
	// If no condition:
	int solutions; // 0,1
	// If solutions=1 and d > 4 (d > 64 for cornacchia4):
	// All solutions are (x,y), (-x,y), (x,-y), (-x,-y).
	cl_I solution_x; // x >= 0
	cl_I solution_y; // y >= 0
	// Constructors.
	cornacchia_t () {}
	cornacchia_t (cl_composite_condition* c) : condition (c) {}
	cornacchia_t (int s) : condition (NULL), solutions (s) {}
	cornacchia_t (int s, const cl_I& x, const cl_I& y) : condition (NULL), solutions (s), solution_x (x), solution_y (y) {}
};
extern const cornacchia_t cornacchia1 (const cl_I& d, const cl_I& p);
extern const cornacchia_t cornacchia4 (const cl_I& d, const cl_I& p);

}  // namespace cln

#endif /* _CL_NUMTHEORY_H */
