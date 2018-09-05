// Univariate Polynomials over the real numbers.

#ifndef _CL_UNIVPOLY_REAL_H
#define _CL_UNIVPOLY_REAL_H

#include "cln/ring.h"
#include "cln/univpoly.h"
#include "cln/number.h"
#include "cln/real_class.h"
#include "cln/integer_class.h"
#include "cln/real_ring.h"

namespace cln {

// Normal univariate polynomials with stricter static typing:
// `cl_R' instead of `cl_ring_element'.

#ifdef notyet

typedef cl_UP_specialized<cl_R> cl_UP_R;
typedef cl_univpoly_specialized_ring<cl_R> cl_univpoly_real_ring;
//typedef cl_heap_univpoly_specialized_ring<cl_R> cl_heap_univpoly_real_ring;

#else

class cl_heap_univpoly_real_ring;

class cl_univpoly_real_ring : public cl_univpoly_ring {
public:
	// Default constructor.
	cl_univpoly_real_ring () : cl_univpoly_ring () {}
	// Copy constructor.
	cl_univpoly_real_ring (const cl_univpoly_real_ring&);
	// Assignment operator.
	cl_univpoly_real_ring& operator= (const cl_univpoly_real_ring&);
	// Automatic dereferencing.
	cl_heap_univpoly_real_ring* operator-> () const
	{ return (cl_heap_univpoly_real_ring*)heappointer; }
};
// Copy constructor and assignment operator.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_univpoly_real_ring,cl_univpoly_ring)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_univpoly_real_ring,cl_univpoly_real_ring)

class cl_UP_R : public cl_UP {
public:
	const cl_univpoly_real_ring& ring () const { return The(cl_univpoly_real_ring)(_ring); }
	// Conversion.
	CL_DEFINE_CONVERTER(cl_ring_element)
	// Destructive modification.
	void set_coeff (uintL index, const cl_R& y);
	void finalize();
	// Evaluation.
	const cl_R operator() (const cl_R& y) const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

class cl_heap_univpoly_real_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
	// High-level operations.
	void fprint (std::ostream& stream, const cl_UP_R& x)
	{
		cl_heap_univpoly_ring::fprint(stream,x);
	}
	bool equal (const cl_UP_R& x, const cl_UP_R& y)
	{
		return cl_heap_univpoly_ring::equal(x,y);
	}
	const cl_UP_R zero ()
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::zero());
	}
	bool zerop (const cl_UP_R& x)
	{
		return cl_heap_univpoly_ring::zerop(x);
	}
	const cl_UP_R plus (const cl_UP_R& x, const cl_UP_R& y)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::plus(x,y));
	}
	const cl_UP_R minus (const cl_UP_R& x, const cl_UP_R& y)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::minus(x,y));
	}
	const cl_UP_R uminus (const cl_UP_R& x)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::uminus(x));
	}
	const cl_UP_R one ()
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::one());
	}
	const cl_UP_R canonhom (const cl_I& x)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::canonhom(x));
	}
	const cl_UP_R mul (const cl_UP_R& x, const cl_UP_R& y)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::mul(x,y));
	}
	const cl_UP_R square (const cl_UP_R& x)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::square(x));
	}
	const cl_UP_R expt_pos (const cl_UP_R& x, const cl_I& y)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::expt_pos(x,y));
	}
	const cl_UP_R scalmul (const cl_R& x, const cl_UP_R& y)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::scalmul(cl_ring_element(cl_R_ring,x),y));
	}
	sintL degree (const cl_UP_R& x)
	{
		return cl_heap_univpoly_ring::degree(x);
	}
	sintL ldegree (const cl_UP_R& x)
	{
		return cl_heap_univpoly_ring::ldegree(x);
	}
	const cl_UP_R monomial (const cl_R& x, uintL e)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::monomial(cl_ring_element(cl_R_ring,x),e));
	}
	const cl_R coeff (const cl_UP_R& x, uintL index)
	{
		return The(cl_R)(cl_heap_univpoly_ring::coeff(x,index));
	}
	const cl_UP_R create (sintL deg)
	{
		return The2(cl_UP_R)(cl_heap_univpoly_ring::create(deg));
	}
	void set_coeff (cl_UP_R& x, uintL index, const cl_R& y)
	{
		cl_heap_univpoly_ring::set_coeff(x,index,cl_ring_element(cl_R_ring,y));
	}
	void finalize (cl_UP_R& x)
	{
		cl_heap_univpoly_ring::finalize(x);
	}
	const cl_R eval (const cl_UP_R& x, const cl_R& y)
	{
		return The(cl_R)(cl_heap_univpoly_ring::eval(x,cl_ring_element(cl_R_ring,y)));
	}
private:
	// No need for any constructors.
	cl_heap_univpoly_real_ring ();
};

// Lookup of polynomial rings.
inline const cl_univpoly_real_ring find_univpoly_ring (const cl_real_ring& r)
{ return The(cl_univpoly_real_ring) (find_univpoly_ring((const cl_ring&)r)); }
inline const cl_univpoly_real_ring find_univpoly_ring (const cl_real_ring& r, const cl_symbol& varname)
{ return The(cl_univpoly_real_ring) (find_univpoly_ring((const cl_ring&)r,varname)); }

// Operations on polynomials.

// Add.
inline const cl_UP_R operator+ (const cl_UP_R& x, const cl_UP_R& y)
	{ return x.ring()->plus(x,y); }

// Negate.
inline const cl_UP_R operator- (const cl_UP_R& x)
	{ return x.ring()->uminus(x); }

// Subtract.
inline const cl_UP_R operator- (const cl_UP_R& x, const cl_UP_R& y)
	{ return x.ring()->minus(x,y); }

// Multiply.
inline const cl_UP_R operator* (const cl_UP_R& x, const cl_UP_R& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
inline const cl_UP_R square (const cl_UP_R& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
inline const cl_UP_R expt_pos (const cl_UP_R& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Scalar multiplication.
#if 0 // less efficient
inline const cl_UP_R operator* (const cl_I& x, const cl_UP_R& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
inline const cl_UP_R operator* (const cl_UP_R& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }
#endif
inline const cl_UP_R operator* (const cl_I& x, const cl_UP_R& y)
	{ return y.ring()->scalmul(x,y); }
inline const cl_UP_R operator* (const cl_UP_R& x, const cl_I& y)
	{ return x.ring()->scalmul(y,x); }
inline const cl_UP_R operator* (const cl_R& x, const cl_UP_R& y)
	{ return y.ring()->scalmul(x,y); }
inline const cl_UP_R operator* (const cl_UP_R& x, const cl_R& y)
	{ return x.ring()->scalmul(y,x); }

// Coefficient.
inline const cl_R coeff (const cl_UP_R& x, uintL index)
	{ return x.ring()->coeff(x,index); }

// Destructive modification.
inline void set_coeff (cl_UP_R& x, uintL index, const cl_R& y)
	{ x.ring()->set_coeff(x,index,y); }
inline void finalize (cl_UP_R& x)
	{ x.ring()->finalize(x); }
inline void cl_UP_R::set_coeff (uintL index, const cl_R& y)
	{ ring()->set_coeff(*this,index,y); }
inline void cl_UP_R::finalize ()
	{ ring()->finalize(*this); }

// Evaluation. (No extension of the base ring allowed here for now.)
inline const cl_R cl_UP_R::operator() (const cl_R& y) const
{
	return ring()->eval(*this,y);
}

// Derivative.
inline const cl_UP_R deriv (const cl_UP_R& x)
	{ return The2(cl_UP_R)(deriv((const cl_UP&)x)); }

#endif

}  // namespace cln

#endif /* _CL_UNIVPOLY_REAL_H */
