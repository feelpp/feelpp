// Univariate Polynomials over the rational numbers.

#ifndef _CL_UNIVPOLY_RATIONAL_H
#define _CL_UNIVPOLY_RATIONAL_H

#include "cln/ring.h"
#include "cln/univpoly.h"
#include "cln/number.h"
#include "cln/rational_class.h"
#include "cln/integer_class.h"
#include "cln/rational_ring.h"

namespace cln {

// Normal univariate polynomials with stricter static typing:
// `cl_RA' instead of `cl_ring_element'.

#ifdef notyet

typedef cl_UP_specialized<cl_RA> cl_UP_RA;
typedef cl_univpoly_specialized_ring<cl_RA> cl_univpoly_rational_ring;
//typedef cl_heap_univpoly_specialized_ring<cl_RA> cl_heap_univpoly_rational_ring;

#else

class cl_heap_univpoly_rational_ring;

class cl_univpoly_rational_ring : public cl_univpoly_ring {
public:
	// Default constructor.
	cl_univpoly_rational_ring () : cl_univpoly_ring () {}
	// Copy constructor.
	cl_univpoly_rational_ring (const cl_univpoly_rational_ring&);
	// Assignment operator.
	cl_univpoly_rational_ring& operator= (const cl_univpoly_rational_ring&);
	// Automatic dereferencing.
	cl_heap_univpoly_rational_ring* operator-> () const
	{ return (cl_heap_univpoly_rational_ring*)heappointer; }
};
// Copy constructor and assignment operator.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_univpoly_rational_ring,cl_univpoly_ring)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_univpoly_rational_ring,cl_univpoly_rational_ring)

class cl_UP_RA : public cl_UP {
public:
	const cl_univpoly_rational_ring& ring () const { return The(cl_univpoly_rational_ring)(_ring); }
	// Conversion.
	CL_DEFINE_CONVERTER(cl_ring_element)
	// Destructive modification.
	void set_coeff (uintL index, const cl_RA& y);
	void finalize();
	// Evaluation.
	const cl_RA operator() (const cl_RA& y) const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

class cl_heap_univpoly_rational_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
	// High-level operations.
	void fprint (std::ostream& stream, const cl_UP_RA& x)
	{
		cl_heap_univpoly_ring::fprint(stream,x);
	}
	bool equal (const cl_UP_RA& x, const cl_UP_RA& y)
	{
		return cl_heap_univpoly_ring::equal(x,y);
	}
	const cl_UP_RA zero ()
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::zero());
	}
	bool zerop (const cl_UP_RA& x)
	{
		return cl_heap_univpoly_ring::zerop(x);
	}
	const cl_UP_RA plus (const cl_UP_RA& x, const cl_UP_RA& y)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::plus(x,y));
	}
	const cl_UP_RA minus (const cl_UP_RA& x, const cl_UP_RA& y)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::minus(x,y));
	}
	const cl_UP_RA uminus (const cl_UP_RA& x)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::uminus(x));
	}
	const cl_UP_RA one ()
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::one());
	}
	const cl_UP_RA canonhom (const cl_I& x)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::canonhom(x));
	}
	const cl_UP_RA mul (const cl_UP_RA& x, const cl_UP_RA& y)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::mul(x,y));
	}
	const cl_UP_RA square (const cl_UP_RA& x)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::square(x));
	}
	const cl_UP_RA expt_pos (const cl_UP_RA& x, const cl_I& y)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::expt_pos(x,y));
	}
	const cl_UP_RA scalmul (const cl_RA& x, const cl_UP_RA& y)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::scalmul(cl_ring_element(cl_RA_ring,x),y));
	}
	sintL degree (const cl_UP_RA& x)
	{
		return cl_heap_univpoly_ring::degree(x);
	}
	sintL ldegree (const cl_UP_RA& x)
	{
		return cl_heap_univpoly_ring::ldegree(x);
	}
	const cl_UP_RA monomial (const cl_RA& x, uintL e)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::monomial(cl_ring_element(cl_RA_ring,x),e));
	}
	const cl_RA coeff (const cl_UP_RA& x, uintL index)
	{
		return The(cl_RA)(cl_heap_univpoly_ring::coeff(x,index));
	}
	const cl_UP_RA create (sintL deg)
	{
		return The2(cl_UP_RA)(cl_heap_univpoly_ring::create(deg));
	}
	void set_coeff (cl_UP_RA& x, uintL index, const cl_RA& y)
	{
		cl_heap_univpoly_ring::set_coeff(x,index,cl_ring_element(cl_RA_ring,y));
	}
	void finalize (cl_UP_RA& x)
	{
		cl_heap_univpoly_ring::finalize(x);
	}
	const cl_RA eval (const cl_UP_RA& x, const cl_RA& y)
	{
		return The(cl_RA)(cl_heap_univpoly_ring::eval(x,cl_ring_element(cl_RA_ring,y)));
	}
private:
	// No need for any constructors.
	cl_heap_univpoly_rational_ring ();
};

// Lookup of polynomial rings.
inline const cl_univpoly_rational_ring find_univpoly_ring (const cl_rational_ring& r)
{ return The(cl_univpoly_rational_ring) (find_univpoly_ring((const cl_ring&)r)); }
inline const cl_univpoly_rational_ring find_univpoly_ring (const cl_rational_ring& r, const cl_symbol& varname)
{ return The(cl_univpoly_rational_ring) (find_univpoly_ring((const cl_ring&)r,varname)); }

// Operations on polynomials.

// Add.
inline const cl_UP_RA operator+ (const cl_UP_RA& x, const cl_UP_RA& y)
	{ return x.ring()->plus(x,y); }

// Negate.
inline const cl_UP_RA operator- (const cl_UP_RA& x)
	{ return x.ring()->uminus(x); }

// Subtract.
inline const cl_UP_RA operator- (const cl_UP_RA& x, const cl_UP_RA& y)
	{ return x.ring()->minus(x,y); }

// Multiply.
inline const cl_UP_RA operator* (const cl_UP_RA& x, const cl_UP_RA& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
inline const cl_UP_RA square (const cl_UP_RA& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
inline const cl_UP_RA expt_pos (const cl_UP_RA& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Scalar multiplication.
#if 0 // less efficient
inline const cl_UP_RA operator* (const cl_I& x, const cl_UP_RA& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
inline const cl_UP_RA operator* (const cl_UP_RA& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }
#endif
inline const cl_UP_RA operator* (const cl_I& x, const cl_UP_RA& y)
	{ return y.ring()->scalmul(x,y); }
inline const cl_UP_RA operator* (const cl_UP_RA& x, const cl_I& y)
	{ return x.ring()->scalmul(y,x); }
inline const cl_UP_RA operator* (const cl_RA& x, const cl_UP_RA& y)
	{ return y.ring()->scalmul(x,y); }
inline const cl_UP_RA operator* (const cl_UP_RA& x, const cl_RA& y)
	{ return x.ring()->scalmul(y,x); }

// Coefficient.
inline const cl_RA coeff (const cl_UP_RA& x, uintL index)
	{ return x.ring()->coeff(x,index); }

// Destructive modification.
inline void set_coeff (cl_UP_RA& x, uintL index, const cl_RA& y)
	{ x.ring()->set_coeff(x,index,y); }
inline void finalize (cl_UP_RA& x)
	{ x.ring()->finalize(x); }
inline void cl_UP_RA::set_coeff (uintL index, const cl_RA& y)
	{ ring()->set_coeff(*this,index,y); }
inline void cl_UP_RA::finalize ()
	{ ring()->finalize(*this); }

// Evaluation. (No extension of the base ring allowed here for now.)
inline const cl_RA cl_UP_RA::operator() (const cl_RA& y) const
{
	return ring()->eval(*this,y);
}

// Derivative.
inline const cl_UP_RA deriv (const cl_UP_RA& x)
	{ return The2(cl_UP_RA)(deriv((const cl_UP&)x)); }

#endif


// Returns the n-th Legendre polynomial (n >= 0).
extern const cl_UP_RA legendre (sintL n);

}  // namespace cln

#endif /* _CL_UNIVPOLY_RATIONAL_H */
