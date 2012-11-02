// Univariate Polynomials over modular integers.

#ifndef _CL_UNIVPOLY_MODINT_H
#define _CL_UNIVPOLY_MODINT_H

#include "cln/ring.h"
#include "cln/univpoly.h"
#include "cln/modinteger.h"
#include "cln/integer_class.h"

namespace cln {

// Normal univariate polynomials with stricter static typing:
// `cl_MI' instead of `cl_ring_element'.

class cl_heap_univpoly_modint_ring;

class cl_univpoly_modint_ring : public cl_univpoly_ring {
public:
	// Default constructor.
	cl_univpoly_modint_ring () : cl_univpoly_ring () {}
	// Copy constructor.
	cl_univpoly_modint_ring (const cl_univpoly_modint_ring&);
	// Assignment operator.
	cl_univpoly_modint_ring& operator= (const cl_univpoly_modint_ring&);
	// Automatic dereferencing.
	cl_heap_univpoly_modint_ring* operator-> () const
	{ return (cl_heap_univpoly_modint_ring*)heappointer; }
};
// Copy constructor and assignment operator.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_univpoly_modint_ring,cl_univpoly_ring)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_univpoly_modint_ring,cl_univpoly_modint_ring)

class cl_UP_MI : public cl_UP {
public:
	const cl_univpoly_modint_ring& ring () const { return The(cl_univpoly_modint_ring)(_ring); }
	// Conversion.
	CL_DEFINE_CONVERTER(cl_ring_element)
	// Destructive modification.
	void set_coeff (uintL index, const cl_MI& y);
	void finalize();
	// Evaluation.
	const cl_MI operator() (const cl_MI& y) const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

class cl_heap_univpoly_modint_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
	const cl_modint_ring& basering () const { return The(cl_modint_ring)(_basering); }
	// High-level operations.
	void fprint (std::ostream& stream, const cl_UP_MI& x)
	{
		cl_heap_univpoly_ring::fprint(stream,x);
	}
	bool equal (const cl_UP_MI& x, const cl_UP_MI& y)
	{
		return cl_heap_univpoly_ring::equal(x,y);
	}
	const cl_UP_MI zero ()
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::zero());
	}
	bool zerop (const cl_UP_MI& x)
	{
		return cl_heap_univpoly_ring::zerop(x);
	}
	const cl_UP_MI plus (const cl_UP_MI& x, const cl_UP_MI& y)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::plus(x,y));
	}
	const cl_UP_MI minus (const cl_UP_MI& x, const cl_UP_MI& y)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::minus(x,y));
	}
	const cl_UP_MI uminus (const cl_UP_MI& x)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::uminus(x));
	}
	const cl_UP_MI one ()
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::one());
	}
	const cl_UP_MI canonhom (const cl_I& x)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::canonhom(x));
	}
	const cl_UP_MI mul (const cl_UP_MI& x, const cl_UP_MI& y)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::mul(x,y));
	}
	const cl_UP_MI square (const cl_UP_MI& x)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::square(x));
	}
	const cl_UP_MI expt_pos (const cl_UP_MI& x, const cl_I& y)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::expt_pos(x,y));
	}
	const cl_UP_MI scalmul (const cl_MI& x, const cl_UP_MI& y)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::scalmul(x,y));
	}
	sintL degree (const cl_UP_MI& x)
	{
		return cl_heap_univpoly_ring::degree(x);
	}
	sintL ldegree (const cl_UP_MI& x)
	{
		return cl_heap_univpoly_ring::ldegree(x);
	}
	const cl_UP_MI monomial (const cl_MI& x, uintL e)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::monomial(x,e));
	}
	const cl_MI coeff (const cl_UP_MI& x, uintL index)
	{
		return The2(cl_MI)(cl_heap_univpoly_ring::coeff(x,index));
	}
	const cl_UP_MI create (sintL deg)
	{
		return The2(cl_UP_MI)(cl_heap_univpoly_ring::create(deg));
	}
	void set_coeff (cl_UP_MI& x, uintL index, const cl_MI& y)
	{
		cl_heap_univpoly_ring::set_coeff(x,index,y);
	}
	void finalize (cl_UP_MI& x)
	{
		cl_heap_univpoly_ring::finalize(x);
	}
	const cl_MI eval (const cl_UP_MI& x, const cl_MI& y)
	{
		return The2(cl_MI)(cl_heap_univpoly_ring::eval(x,y));
	}
private:
	// No need for any constructors.
	cl_heap_univpoly_modint_ring ();
};

// Lookup of polynomial rings.
inline const cl_univpoly_modint_ring find_univpoly_ring (const cl_modint_ring& r)
{ return The(cl_univpoly_modint_ring) (find_univpoly_ring((const cl_ring&)r)); }
inline const cl_univpoly_modint_ring find_univpoly_ring (const cl_modint_ring& r, const cl_symbol& varname)
{ return The(cl_univpoly_modint_ring) (find_univpoly_ring((const cl_ring&)r,varname)); }

// Operations on polynomials.

// Add.
inline const cl_UP_MI operator+ (const cl_UP_MI& x, const cl_UP_MI& y)
	{ return x.ring()->plus(x,y); }

// Negate.
inline const cl_UP_MI operator- (const cl_UP_MI& x)
	{ return x.ring()->uminus(x); }

// Subtract.
inline const cl_UP_MI operator- (const cl_UP_MI& x, const cl_UP_MI& y)
	{ return x.ring()->minus(x,y); }

// Multiply.
inline const cl_UP_MI operator* (const cl_UP_MI& x, const cl_UP_MI& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
inline const cl_UP_MI square (const cl_UP_MI& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
inline const cl_UP_MI expt_pos (const cl_UP_MI& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Scalar multiplication.
#if 0 // less efficient
inline const cl_UP_MI operator* (const cl_I& x, const cl_UP_MI& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
inline const cl_UP_MI operator* (const cl_UP_MI& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }
#endif
inline const cl_UP_MI operator* (const cl_I& x, const cl_UP_MI& y)
	{ return y.ring()->scalmul(y.ring()->basering()->canonhom(x),y); }
inline const cl_UP_MI operator* (const cl_UP_MI& x, const cl_I& y)
	{ return x.ring()->scalmul(x.ring()->basering()->canonhom(y),x); }
inline const cl_UP_MI operator* (const cl_MI& x, const cl_UP_MI& y)
	{ return y.ring()->scalmul(x,y); }
inline const cl_UP_MI operator* (const cl_UP_MI& x, const cl_MI& y)
	{ return x.ring()->scalmul(y,x); }

// Coefficient.
inline const cl_MI coeff (const cl_UP_MI& x, uintL index)
	{ return x.ring()->coeff(x,index); }

// Destructive modification.
inline void set_coeff (cl_UP_MI& x, uintL index, const cl_MI& y)
	{ x.ring()->set_coeff(x,index,y); }
inline void finalize (cl_UP_MI& x)
	{ x.ring()->finalize(x); }
inline void cl_UP_MI::set_coeff (uintL index, const cl_MI& y)
	{ ring()->set_coeff(*this,index,y); }
inline void cl_UP_MI::finalize ()
	{ ring()->finalize(*this); }

// Evaluation. (No extension of the base ring allowed here for now.)
inline const cl_MI cl_UP_MI::operator() (const cl_MI& y) const
{
	return ring()->eval(*this,y);
}

// Derivative.
inline const cl_UP_MI deriv (const cl_UP_MI& x)
	{ return The2(cl_UP_MI)(deriv((const cl_UP&)x)); }

}  // namespace cln

#endif /* _CL_UNIVPOLY_MODINT_H */
