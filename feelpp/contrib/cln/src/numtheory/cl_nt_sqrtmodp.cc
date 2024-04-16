// sqrt_mod_p().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/numtheory.h"


// Implementation.

#include "integer/cl_I.h"
#include "cln/exception.h"

#undef floor
#include <cmath>
#define floor cln_floor

// MacOS X does "#define _R 0x00040000L".  Grr...
#undef _R

namespace cln {

// Algorithm 1 (for very small p only):
// Try different values.
// Assume p is prime and a nonzero square in Z/pZ.
static uint32 search_sqrt (uint32 p, uint32 a)
{
	var uint32 x = 1;
	var uint32 x2 = 1;
	loop {
		// 0 < x <= p/2, x2 = x^2 mod p.
		if (x2 == a)
			return x;
		x2 += x; x++; x2 += x;
		if (x2 >= p)
			x2 -= p;
	}
}

// Algorithm 2 (for p > 2 only):
// Cantor-Zassenhaus.
// [Beth et al.: Computer Algebra, 1988, Kapitel 5.3.3.]
// [Cohen, A Course in Computational Algebraic Number Theory,
//  Section 3.4.4., Algorithm 3.4.6.]
// Input: R = Z/pZ with p>2, and a (nonzero square in R).
static const sqrt_mod_p_t cantor_zassenhaus_sqrt (const cl_modint_ring& R, const cl_MI& a);
	// Compute in the polynomial ring R[X]/(X^2-a).
	struct pol2 {
		// A polynomial c0+c1*X mod (X^2-a)
		cl_MI c0;
		cl_MI c1;
		// Constructor.
		pol2 (const cl_MI& _c0, const cl_MI& _c1) : c0 (_c0), c1 (_c1) {}
	};
	struct pol2ring {
		const cl_modint_ring& R;
		const cl_MI& a;
		const pol2 zero ()
		{
			return pol2(R->zero(),R->zero());
		}
		const pol2 one ()
		{
			return pol2(R->one(),R->zero());
		}
		const pol2 plus (const pol2& u, const pol2& v)
		{
			return pol2(u.c0+v.c0, u.c1+v.c1);
		}
		const pol2 minus (const pol2& u, const pol2& v)
		{
			return pol2(u.c0-v.c0, u.c1-v.c1);
		}
		const pol2 mul (const pol2& u, const pol2& v)
		{
			return pol2(u.c0*v.c0+u.c1*v.c1*a, u.c0*v.c1+u.c1*v.c0);
		}
		const pol2 square (const pol2& u)
		{
			return pol2(cln::square(u.c0) + cln::square(u.c1)*a, (u.c0*u.c1)<<1);
		}
		const pol2 expt_pos (const pol2& x, const cl_I& y)
		{
			// Right-Left Binary, [Cohen, Algorithm 1.2.1.]
			var pol2 a = x;
			var cl_I b = y;
			while (!oddp(b)) { a = square(a); b = b = b >> 1; } // a^b = x^y
			var pol2 c = a;
			until (eq(b,1)) {
				b = b >> 1;
				a = square(a);
				// a^b*c = x^y
				if (oddp(b))
					c = mul(a,c);
			}
			return c;
		}
		const pol2 random ()
		{
			return pol2(R->random(),R->random());
		}
		// Computes the degree of gcd(u(X),X^2-a) and, if it is 1,
		// also the zero if this polynomial of degree 1.
		struct gcd_result {
			cl_composite_condition* condition;
			int gcd_degree;
			cl_MI solution;
			// Constructors.
			gcd_result (cl_composite_condition* c) : condition (c) {}
			gcd_result (int deg) : condition (NULL), gcd_degree (deg) {}
			gcd_result (int deg, const cl_MI& sol) : condition (NULL), gcd_degree (deg), solution (sol) {}
		};
		const gcd_result gcd (const pol2& u)
		{
			if (zerop(u.c1)) {
				// constant polynomial u(X)
				if (zerop(u.c0))
					return gcd_result(2);
				else
					return gcd_result(0);
			}
			// u(X) = c0 + c1*X has zero -c0/c1.
			var cl_MI_x c1inv = R->recip(u.c1);
			if (c1inv.condition)
				return c1inv.condition;
			var cl_MI z = -u.c0*c1inv;
			if (cln::square(z) == a)
				return gcd_result(1,z);
			else
				return gcd_result(0);
		}
		// Constructor.
		pol2ring (const cl_modint_ring& _R, const cl_MI& _a) : R (_R), a (_a) {}
	};
static const sqrt_mod_p_t cantor_zassenhaus_sqrt (const cl_modint_ring& R, const cl_MI& a)
{
	var pol2ring PR = pol2ring(R,a);
	var cl_I& p = R->modulus;
	// Assuming p is a prime, then R[X]/(X^2-a) is the direct product of
	// two rings R[X]/(X-sqrt(a)), each being isomorphic to R. Thus taking
	// a (p-1)/2-th power in this ring will return one of (0,+1,-1) in
	// each ring, with independent probabilities (1/p, (p-1)/2p, (p-1)/2p).
	// For any polynomial u(X), setting v(X) := u(X)^((p-1)/2) yields
	// gcd(u(X),X^2-a) * gcd(v(X)-1,X^2-a) * gcd(v(X)+1,X^2-a) = X^2-a.
	// If p is not prime, all of these gcd's are likely to be 1.
	var cl_I e = (p-1) >> 1;
	loop {
		// Choose a random polynomial u(X) in the ring.
		var pol2 u = PR.random();
		// Compute v(X) = u(X)^((p-1)/2).
		var pol2 v = PR.expt_pos(u,e);
		// Compute the three gcds.
		var pol2ring::gcd_result g1 = PR.gcd(PR.minus(v,PR.one()));
		if (g1.condition)
			return g1.condition;
		if (g1.gcd_degree == 1)
			return sqrt_mod_p_t(2,g1.solution,-g1.solution);
		if (g1.gcd_degree == 2)
			continue;
		var pol2ring::gcd_result g2 = PR.gcd(PR.plus(v,PR.one()));
		if (g2.condition)
			return g2.condition;
		if (g2.gcd_degree == 1)
			return sqrt_mod_p_t(2,g2.solution,-g2.solution);
		if (g2.gcd_degree == 2)
			continue;
		var pol2ring::gcd_result g3 = PR.gcd(u);
		if (g3.condition)
			return g3.condition;
		if (g3.gcd_degree == 1)
			return sqrt_mod_p_t(2,g3.solution,-g3.solution);
		if (g1.gcd_degree + g2.gcd_degree + g3.gcd_degree < 2)
			// If the sum of the degrees of the gcd is != 2,
			// p cannot be prime.
			return new cl_composite_condition(p);
	}
}

// Algorithm 3 (for p > 2 only):
// Tonelli-Shanks.
// [Cohen, A Course in Computational Algebraic Number Theory,
//  Section 1.5.1., Algorithm 1.5.1.]
static const sqrt_mod_p_t tonelli_shanks_sqrt (const cl_modint_ring& R, const cl_MI& a)
{
	// Idea:
	// Write p-1 = 2^e*m, m odd. G = (Z/pZ)^* (cyclic of order p-1) has
	// subgroups G_0 < G_1 < ... < G_e, G_j of order 2^j. (G_e is called
	// the "2-Sylow subgroup" of G.) More precisely
	//          G_j = { x in (Z/pZ)^* : x^(2^j) = 1 },
	//        G/G_j = { x^(2^j) : x in (Z/pZ)^* }.
	// We compute the square root of a first in G/G_e, then lift it to
	// G/G_(e-1), etc., up to G/G_0.
	// Start with b = a^((m+1)/2), then (a^-1*b^2)^(2^e) = 1, i.e.
	// a = b^2 in G/G_e.
	// Lifting from G/G_j to G/G_(j-1) is easy: Assume a = b^2 in G/G_j.
	// If a = b^2 in G/G_(j-1), then nothing needs to be done. Else
	// a^-1*b^2 is in G_j \ G_(j-1). If j=e, a^-1*b^2 is a non-square
	// mod p, hence a is a non-square as well, contradiction. If j<e,
	// take h in G_(j+1) \ G_j, so that h^2 in G_j \ G_(j-1), and
	// a^-1*b^2*h^2 is in G_(j-1). So multiply b with h.
	var cl_I& p = R->modulus;
	var uintC e = ord2(p-1);
	var cl_I m = (p-1) >> e;
	// p-1 = 2^e*m, m odd.
	// We will have the invariant c = a^-1*b^2 in G/G_j.
	var uintC j = e;
	// Initialize b = a^((m+1)/2), c = a^m, but avoid to divide by a.
	var cl_MI c = R->expt_pos(a,(m-1)>>1);
	var cl_MI b = R->mul(a,c);
	c = R->mul(b,c);
	// Find h in G_e \ G_(e-1): h = h'^m, where h' is any non-square.
	var cl_MI h;
	if (e==1)
		h = - R->one();
	else {
		// Since this computation is a bit costly, we cache its result
		// on the ring's property list.
		static const cl_symbol key = (cl_symbol)(cl_string)"generator of 2-Sylow subgroup of (Z/pZ)^*";
		struct cl_sylow2gen_property : public cl_property {
			SUBCLASS_cl_property();
		public:
			cl_I h_rep;
			// Constructor.
			cl_sylow2gen_property (const cl_symbol& k, const cl_MI& h) : cl_property (k), h_rep (h.rep) {}
		};
		var cl_sylow2gen_property* prop = (cl_sylow2gen_property*) R->get_property(key);
		if (prop)
			h = cl_MI(R,prop->h_rep);
		else {
			do { h = R->random(); }
			   until (jacobi(R->retract(h),p) == -1);
			h = R->expt_pos(h,m);
			R->add_property(new cl_sylow2gen_property(key,h));
		}
	}
	do {
		// Now c = a^-1*b^2 in G_j, h in G_j \ G_(j-1).
		// Determine the smallest i such that c in G_i.
		var uintC i = 0;
		var cl_MI ci = c; // c_i = c^(2^i)
		for ( ; i < j; i++, ci = R->square(ci))
			if (ci == R->one())
				break;
		if (i==j)
			// Some problem: if j=e, a non-square, if j<e, the
			// previous iteration didn't do its job correctly.
			// Indicates that p is not prime.
			return new cl_composite_condition(p);
		// OK, i < j.
		for (var uintC count = j-i-1; count > 0; count--)
			h = R->square(h);
		// Now h in G_(i+1) \ G_i.
		b = R->mul(b,h);
		h = R->square(h);
		c = R->mul(c,h);
		// Now c = a^-1*b^2 in G_(i-1), h in G_i \ G_(i-1).
		j = i;
	} while (j > 0);
	if (R->square(b) != a)
		// Problem again.
		return new cl_composite_condition(p);
	return sqrt_mod_p_t(2,b,-b);
}

// Break-Even-Points (on a i486 with 33 MHz):
// Algorithm 1 fastest for p < 1500..2000
// Algorithm 3 generally fastest for p > 2000.
// But the running time of algorithm 3 is proportional to e^2.
// For large e, algorithm 2 becomes faster.
// l=50 bits: for e >= 40
// l=100 bits: for e >= 55
// l=200 bits: for e >= 80
// l=400 bits: for e >= 130
// in general something like  e > l/(log(l)/(2*log(2))-1).

const sqrt_mod_p_t sqrt_mod_p (const cl_modint_ring& R, const cl_MI& a)
{
	if (!(a.ring() == R)) throw runtime_exception();
	var cl_I& p = R->modulus;
	var cl_I aa = R->retract(a);
	switch (jacobi(aa,p)) {
		case -1: // no solution
			return sqrt_mod_p_t(0);
		case 0: // gcd(aa,p) > 1
			if (zerop(a))
				// one solution
				return sqrt_mod_p_t(1,a);
			else
				// found factor of p
				return new cl_composite_condition(p,gcd(aa,p));
		case 1: // two solutions
			break;
	}
	if (p < 2000) {
		// Algorithm 1.
		var cl_I x1 = search_sqrt(cl_I_to_UL(p),cl_I_to_UL(aa));
		var cl_I x2 = p-x1;
		if (x1==x2) // can only happen when p = 2
			return sqrt_mod_p_t(1,R->canonhom(x1));
		else
			return sqrt_mod_p_t(2,R->canonhom(x1),R->canonhom(x2));
	}
	var uintC l = integer_length(p);
	var uintC e = ord2(p-1);
	//if (e > 30 && e > l/(::log((double)l)*0.72-1))
	if (e > 30 && e > l/(::log((double)l)*0.92-2.41))
		// Algorithm 2.
		return cantor_zassenhaus_sqrt(R,a);
	else
		// Algorithm 3.
		return tonelli_shanks_sqrt(R,a);
}

}  // namespace cln
