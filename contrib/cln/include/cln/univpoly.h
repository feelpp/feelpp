// Univariate Polynomials.

#ifndef _CL_UNIVPOLY_H
#define _CL_UNIVPOLY_H

#include "cln/object.h"
#include "cln/ring.h"
#include "cln/malloc.h"
#include "cln/proplist.h"
#include "cln/symbol.h"
#include "cln/V.h"
#include "cln/io.h"

namespace cln {

// To protect against mixing elements of different polynomial rings, every
// polynomial carries its ring in itself.

class cl_heap_univpoly_ring;

class cl_univpoly_ring : public cl_ring {
public:
	// Default constructor.
	cl_univpoly_ring ();
	// Constructor. Takes a cl_heap_univpoly_ring*, increments its refcount.
	cl_univpoly_ring (cl_heap_univpoly_ring* r);
	// Private constructor. Doesn't increment the refcount.
	cl_univpoly_ring (cl_private_thing);
	// Copy constructor.
	cl_univpoly_ring (const cl_univpoly_ring&);
	// Assignment operator.
	cl_univpoly_ring& operator= (const cl_univpoly_ring&);
	// Automatic dereferencing.
	cl_heap_univpoly_ring* operator-> () const
	{ return (cl_heap_univpoly_ring*)heappointer; }
};
// Copy constructor and assignment operator.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_univpoly_ring,cl_ring)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_univpoly_ring,cl_univpoly_ring)

// Normal constructor for `cl_univpoly_ring'.
inline cl_univpoly_ring::cl_univpoly_ring (cl_heap_univpoly_ring* r)
	: cl_ring ((cl_private_thing) (cl_inc_pointer_refcount((cl_heap*)r), r)) {}
// Private constructor for `cl_univpoly_ring'.
inline cl_univpoly_ring::cl_univpoly_ring (cl_private_thing p)
	: cl_ring (p) {}

// Operations on univariate polynomial rings.

inline bool operator== (const cl_univpoly_ring& R1, const cl_univpoly_ring& R2)
{ return (R1.pointer == R2.pointer); }
inline bool operator!= (const cl_univpoly_ring& R1, const cl_univpoly_ring& R2)
{ return (R1.pointer != R2.pointer); }
inline bool operator== (const cl_univpoly_ring& R1, cl_heap_univpoly_ring* R2)
{ return (R1.pointer == R2); }
inline bool operator!= (const cl_univpoly_ring& R1, cl_heap_univpoly_ring* R2)
{ return (R1.pointer != R2); }

// Representation of a univariate polynomial.

class _cl_UP /* cf. _cl_ring_element */ {
public:
	cl_gcpointer rep; // vector of coefficients, a cl_V_any
	// Default constructor.
	_cl_UP ();
public: /* ugh */
	// Constructor.
	_cl_UP (const cl_heap_univpoly_ring* R, const cl_V_any& r) : rep (as_cl_private_thing(r)) { (void)R; }
	_cl_UP (const cl_univpoly_ring& R, const cl_V_any& r) : rep (as_cl_private_thing(r)) { (void)R; }
public:
	// Conversion.
	CL_DEFINE_CONVERTER(_cl_ring_element)
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

class cl_UP /* cf. cl_ring_element */ : public _cl_UP {
protected:
	cl_univpoly_ring _ring;	// polynomial ring (references the base ring)
public:
	const cl_univpoly_ring& ring () const { return _ring; }
private:
	// Default constructor.
	cl_UP ();
public: /* ugh */
	// Constructor.
	cl_UP (const cl_univpoly_ring& R, const cl_V_any& r)
		: _cl_UP (R,r), _ring (R) {}
	cl_UP (const cl_univpoly_ring& R, const _cl_UP& r)
		: _cl_UP (r), _ring (R) {}
public:
	// Conversion.
	CL_DEFINE_CONVERTER(cl_ring_element)
	// Destructive modification.
	void set_coeff (uintL index, const cl_ring_element& y);
	void finalize();
	// Evaluation.
	const cl_ring_element operator() (const cl_ring_element& y) const;
	// Debugging output.
	void debug_print () const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};


// Ring operations.

struct _cl_univpoly_setops /* cf. _cl_ring_setops */ {
	// print
	void (* fprint) (cl_heap_univpoly_ring* R, std::ostream& stream, const _cl_UP& x);
	// equality
	// (Be careful: This is not well-defined for polynomials with
	// floating-point coefficients.)
	bool (* equal) (cl_heap_univpoly_ring* R, const _cl_UP& x, const _cl_UP& y);
};
struct _cl_univpoly_addops /* cf. _cl_ring_addops */ {
	// 0
	const _cl_UP (* zero) (cl_heap_univpoly_ring* R);
	bool (* zerop) (cl_heap_univpoly_ring* R, const _cl_UP& x);
	// x+y
	const _cl_UP (* plus) (cl_heap_univpoly_ring* R, const _cl_UP& x, const _cl_UP& y);
	// x-y
	const _cl_UP (* minus) (cl_heap_univpoly_ring* R, const _cl_UP& x, const _cl_UP& y);
	// -x
	const _cl_UP (* uminus) (cl_heap_univpoly_ring* R, const _cl_UP& x);
};
struct _cl_univpoly_mulops /* cf. _cl_ring_mulops */ {
	// 1
	const _cl_UP (* one) (cl_heap_univpoly_ring* R);
	// canonical homomorphism
	const _cl_UP (* canonhom) (cl_heap_univpoly_ring* R, const cl_I& x);
	// x*y
	const _cl_UP (* mul) (cl_heap_univpoly_ring* R, const _cl_UP& x, const _cl_UP& y);
	// x^2
	const _cl_UP (* square) (cl_heap_univpoly_ring* R, const _cl_UP& x);
	// x^y, y Integer >0
	const _cl_UP (* expt_pos) (cl_heap_univpoly_ring* R, const _cl_UP& x, const cl_I& y);
};
struct _cl_univpoly_modulops {
	// scalar multiplication x*y
	const _cl_UP (* scalmul) (cl_heap_univpoly_ring* R, const cl_ring_element& x, const _cl_UP& y);
};
struct _cl_univpoly_polyops {
	// degree
	sintL (* degree) (cl_heap_univpoly_ring* R, const _cl_UP& x);
	// low degree
	sintL (* ldegree) (cl_heap_univpoly_ring* R, const _cl_UP& x);
	// monomial
	const _cl_UP (* monomial) (cl_heap_univpoly_ring* R, const cl_ring_element& x, uintL e);
	// coefficient (0 if index>degree)
	const cl_ring_element (* coeff) (cl_heap_univpoly_ring* R, const _cl_UP& x, uintL index);
	// create new polynomial, bounded degree
	const _cl_UP (* create) (cl_heap_univpoly_ring* R, sintL deg);
	// set coefficient in new polynomial
	void (* set_coeff) (cl_heap_univpoly_ring* R, _cl_UP& x, uintL index, const cl_ring_element& y);
	// finalize polynomial
	void (* finalize) (cl_heap_univpoly_ring* R, _cl_UP& x);
	// evaluate, substitute an element of R
	const cl_ring_element (* eval) (cl_heap_univpoly_ring* R, const _cl_UP& x, const cl_ring_element& y);
};
  typedef const _cl_univpoly_setops  cl_univpoly_setops;
  typedef const _cl_univpoly_addops  cl_univpoly_addops;
  typedef const _cl_univpoly_mulops  cl_univpoly_mulops;
  typedef const _cl_univpoly_modulops  cl_univpoly_modulops;
  typedef const _cl_univpoly_polyops  cl_univpoly_polyops;

// Representation of a univariate polynomial ring.

class cl_heap_univpoly_ring /* cf. cl_heap_ring */ : public cl_heap {
	SUBCLASS_cl_heap_ring()
private:
	cl_property_list properties;
protected:
	cl_univpoly_setops* setops;
	cl_univpoly_addops* addops;
	cl_univpoly_mulops* mulops;
	cl_univpoly_modulops* modulops;
	cl_univpoly_polyops* polyops;
protected:
	cl_ring _basering;	// the coefficients are elements of this ring
public:
	const cl_ring& basering () const { return _basering; }
public:
	// Low-level operations.
	void _fprint (std::ostream& stream, const _cl_UP& x)
		{ setops->fprint(this,stream,x); }
	bool _equal (const _cl_UP& x, const _cl_UP& y)
		{ return setops->equal(this,x,y); }
	const _cl_UP _zero ()
		{ return addops->zero(this); }
	bool _zerop (const _cl_UP& x)
		{ return addops->zerop(this,x); }
	const _cl_UP _plus (const _cl_UP& x, const _cl_UP& y)
		{ return addops->plus(this,x,y); }
	const _cl_UP _minus (const _cl_UP& x, const _cl_UP& y)
		{ return addops->minus(this,x,y); }
	const _cl_UP _uminus (const _cl_UP& x)
		{ return addops->uminus(this,x); }
	const _cl_UP _one ()
		{ return mulops->one(this); }
	const _cl_UP _canonhom (const cl_I& x)
		{ return mulops->canonhom(this,x); }
	const _cl_UP _mul (const _cl_UP& x, const _cl_UP& y)
		{ return mulops->mul(this,x,y); }
	const _cl_UP _square (const _cl_UP& x)
		{ return mulops->square(this,x); }
	const _cl_UP _expt_pos (const _cl_UP& x, const cl_I& y)
		{ return mulops->expt_pos(this,x,y); }
	const _cl_UP _scalmul (const cl_ring_element& x, const _cl_UP& y)
		{ return modulops->scalmul(this,x,y); }
	sintL _degree (const _cl_UP& x)
		{ return polyops->degree(this,x); }
	sintL _ldegree (const _cl_UP& x)
		{ return polyops->ldegree(this,x); }
	const _cl_UP _monomial (const cl_ring_element& x, uintL e)
		{ return polyops->monomial(this,x,e); }
	const cl_ring_element _coeff (const _cl_UP& x, uintL index)
		{ return polyops->coeff(this,x,index); }
	const _cl_UP _create (sintL deg)
		{ return polyops->create(this,deg); }
	void _set_coeff (_cl_UP& x, uintL index, const cl_ring_element& y)
		{ polyops->set_coeff(this,x,index,y); }
	void _finalize (_cl_UP& x)
		{ polyops->finalize(this,x); }
	const cl_ring_element _eval (const _cl_UP& x, const cl_ring_element& y)
		{ return polyops->eval(this,x,y); }
	// High-level operations.
	void fprint (std::ostream& stream, const cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		_fprint(stream,x);
	}
	bool equal (const cl_UP& x, const cl_UP& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return _equal(x,y);
	}
	const cl_UP zero ()
	{
		return cl_UP(this,_zero());
	}
	bool zerop (const cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _zerop(x);
	}
	const cl_UP plus (const cl_UP& x, const cl_UP& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_UP(this,_plus(x,y));
	}
	const cl_UP minus (const cl_UP& x, const cl_UP& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_UP(this,_minus(x,y));
	}
	const cl_UP uminus (const cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_UP(this,_uminus(x));
	}
	const cl_UP one ()
	{
		return cl_UP(this,_one());
	}
	const cl_UP canonhom (const cl_I& x)
	{
		return cl_UP(this,_canonhom(x));
	}
	const cl_UP mul (const cl_UP& x, const cl_UP& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_UP(this,_mul(x,y));
	}
	const cl_UP square (const cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_UP(this,_square(x));
	}
	const cl_UP expt_pos (const cl_UP& x, const cl_I& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_UP(this,_expt_pos(x,y));
	}
	const cl_UP scalmul (const cl_ring_element& x, const cl_UP& y)
	{
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_UP(this,_scalmul(x,y));
	}
	sintL degree (const cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _degree(x);
	}
	sintL ldegree (const cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _ldegree(x);
	}
	const cl_UP monomial (const cl_ring_element& x, uintL e)
	{
		return cl_UP(this,_monomial(x,e));
	}
	const cl_ring_element coeff (const cl_UP& x, uintL index)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _coeff(x,index);
	}
	const cl_UP create (sintL deg)
	{
		return cl_UP(this,_create(deg));
	}
	void set_coeff (cl_UP& x, uintL index, const cl_ring_element& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		_set_coeff(x,index,y);
	}
	void finalize (cl_UP& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		_finalize(x);
	}
	const cl_ring_element eval (const cl_UP& x, const cl_ring_element& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _eval(x,y);
	}
	// Property operations.
	cl_property* get_property (const cl_symbol& key)
		{ return properties.get_property(key); }
	void add_property (cl_property* new_property)
		{ properties.add_property(new_property); }
// Constructor.
	cl_heap_univpoly_ring (const cl_ring& r, cl_univpoly_setops*, cl_univpoly_addops*, cl_univpoly_mulops*, cl_univpoly_modulops*, cl_univpoly_polyops*);
	~cl_heap_univpoly_ring () {}
};
#define SUBCLASS_cl_heap_univpoly_ring() \
  SUBCLASS_cl_heap_ring()


// Lookup or create the "standard" univariate polynomial ring over a ring r.
extern const cl_univpoly_ring find_univpoly_ring (const cl_ring& r);

// Lookup or create a univariate polynomial ring with a named variable over r.
extern const cl_univpoly_ring find_univpoly_ring (const cl_ring& r, const cl_symbol& varname);

class cl_UP_init_helper
{
	static int count;
public:
	cl_UP_init_helper();
	~cl_UP_init_helper();
};
static cl_UP_init_helper cl_UP_init_helper_instance;


// Operations on polynomials.

// Output.
inline void fprint (std::ostream& stream, const cl_UP& x)
	{ x.ring()->fprint(stream,x); }
CL_DEFINE_PRINT_OPERATOR(cl_UP)

// Add.
inline const cl_UP operator+ (const cl_UP& x, const cl_UP& y)
	{ return x.ring()->plus(x,y); }

// Negate.
inline const cl_UP operator- (const cl_UP& x)
	{ return x.ring()->uminus(x); }

// Subtract.
inline const cl_UP operator- (const cl_UP& x, const cl_UP& y)
	{ return x.ring()->minus(x,y); }

// Equality.
inline bool operator== (const cl_UP& x, const cl_UP& y)
	{ return x.ring()->equal(x,y); }
inline bool operator!= (const cl_UP& x, const cl_UP& y)
	{ return !x.ring()->equal(x,y); }

// Compare against 0.
inline bool zerop (const cl_UP& x)
	{ return x.ring()->zerop(x); }

// Multiply.
inline const cl_UP operator* (const cl_UP& x, const cl_UP& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
inline const cl_UP square (const cl_UP& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
inline const cl_UP expt_pos (const cl_UP& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Scalar multiplication.
#if 0 // less efficient
inline const cl_UP operator* (const cl_I& x, const cl_UP& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
inline const cl_UP operator* (const cl_UP& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }
#endif
inline const cl_UP operator* (const cl_I& x, const cl_UP& y)
	{ return y.ring()->scalmul(y.ring()->basering()->canonhom(x),y); }
inline const cl_UP operator* (const cl_UP& x, const cl_I& y)
	{ return x.ring()->scalmul(x.ring()->basering()->canonhom(y),x); }
inline const cl_UP operator* (const cl_ring_element& x, const cl_UP& y)
	{ return y.ring()->scalmul(x,y); }
inline const cl_UP operator* (const cl_UP& x, const cl_ring_element& y)
	{ return x.ring()->scalmul(y,x); }

// Degree.
inline sintL degree (const cl_UP& x)
	{ return x.ring()->degree(x); }

// Low degree.
inline sintL ldegree (const cl_UP& x)
	{ return x.ring()->ldegree(x); }

// Coefficient.
inline const cl_ring_element coeff (const cl_UP& x, uintL index)
	{ return x.ring()->coeff(x,index); }

// Destructive modification.
inline void set_coeff (cl_UP& x, uintL index, const cl_ring_element& y)
	{ x.ring()->set_coeff(x,index,y); }
inline void finalize (cl_UP& x)
	{ x.ring()->finalize(x); }
inline void cl_UP::set_coeff (uintL index, const cl_ring_element& y)
	{ ring()->set_coeff(*this,index,y); }
inline void cl_UP::finalize ()
	{ ring()->finalize(*this); }

// Evaluation. (No extension of the base ring allowed here for now.)
inline const cl_ring_element cl_UP::operator() (const cl_ring_element& y) const
{
	return ring()->eval(*this,y);
}

// Derivative.
extern const cl_UP deriv (const cl_UP& x);


// Ring of uninitialized elements.
// Any operation results in a run-time error.

extern const cl_univpoly_ring cl_no_univpoly_ring;
extern cl_class cl_class_no_univpoly_ring;

class cl_UP_no_ring_init_helper
{
	static int count;
public:
	cl_UP_no_ring_init_helper();
	~cl_UP_no_ring_init_helper();
};
static cl_UP_no_ring_init_helper cl_UP_no_ring_init_helper_instance;

inline cl_univpoly_ring::cl_univpoly_ring ()
	: cl_ring (as_cl_private_thing(cl_no_univpoly_ring)) {}
inline _cl_UP::_cl_UP ()
	: rep ((cl_private_thing) cl_combine(cl_FN_tag,0)) {}
inline cl_UP::cl_UP ()
	: _cl_UP (), _ring () {}


// Debugging support.
#ifdef CL_DEBUG
extern int cl_UP_debug_module;
CL_FORCE_LINK(cl_UP_debug_dummy, cl_UP_debug_module)
#endif

}  // namespace cln

#endif /* _CL_UNIVPOLY_H */

namespace cln {

// Templates for univariate polynomials of complex/real/rational/integers.

#ifdef notyet
// Unfortunately, this is not usable now, because of gcc-2.7 bugs:
// - A template inline function is not inline in the first function that
//   uses it.
// - Argument matching bug: User-defined conversions are not tried (or
//   tried with too low priority) for template functions w.r.t. normal
//   functions. For example, a call expt_pos(cl_UP_specialized<cl_N>,int)
//   is compiled as expt_pos(const cl_UP&, const cl_I&) instead of
//   expt_pos(const cl_UP_specialized<cl_N>&, const cl_I&).
// It will, however, be usable when gcc-2.8 is released.

#if defined(_CL_UNIVPOLY_COMPLEX_H) || defined(_CL_UNIVPOLY_REAL_H) || defined(_CL_UNIVPOLY_RATIONAL_H) || defined(_CL_UNIVPOLY_INTEGER_H)
#ifndef _CL_UNIVPOLY_AUX_H

// Normal univariate polynomials with stricter static typing:
// `class T' instead of `cl_ring_element'.

template <class T> class cl_univpoly_specialized_ring;
template <class T> class cl_UP_specialized;
template <class T> class cl_heap_univpoly_specialized_ring;

template <class T>
class cl_univpoly_specialized_ring : public cl_univpoly_ring {
public:
	// Default constructor.
	cl_univpoly_specialized_ring () : cl_univpoly_ring () {}
	// Copy constructor.
	cl_univpoly_specialized_ring (const cl_univpoly_specialized_ring&);
	// Assignment operator.
	cl_univpoly_specialized_ring& operator= (const cl_univpoly_specialized_ring&);
	// Automatic dereferencing.
	cl_heap_univpoly_specialized_ring<T>* operator-> () const
	{ return (cl_heap_univpoly_specialized_ring<T>*)heappointer; }
};
// Copy constructor and assignment operator.
template <class T>
_CL_DEFINE_COPY_CONSTRUCTOR2(cl_univpoly_specialized_ring<T>,cl_univpoly_specialized_ring,cl_univpoly_ring)
template <class T>
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_univpoly_specialized_ring<T>,cl_univpoly_specialized_ring<T>)

template <class T>
class cl_UP_specialized : public cl_UP {
public:
	const cl_univpoly_specialized_ring<T>& ring () const { return The(cl_univpoly_specialized_ring<T>)(_ring); }
	// Conversion.
	CL_DEFINE_CONVERTER(cl_ring_element)
	// Destructive modification.
	void set_coeff (uintL index, const T& y);
	void finalize();
	// Evaluation.
	const T operator() (const T& y) const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

template <class T>
class cl_heap_univpoly_specialized_ring : public cl_heap_univpoly_ring {
	SUBCLASS_cl_heap_univpoly_ring()
	// High-level operations.
	void fprint (std::ostream& stream, const cl_UP_specialized<T>& x)
	{
		cl_heap_univpoly_ring::fprint(stream,x);
	}
	bool equal (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{
		return cl_heap_univpoly_ring::equal(x,y);
	}
	const cl_UP_specialized<T> zero ()
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::zero());
	}
	bool zerop (const cl_UP_specialized<T>& x)
	{
		return cl_heap_univpoly_ring::zerop(x);
	}
	const cl_UP_specialized<T> plus (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::plus(x,y));
	}
	const cl_UP_specialized<T> minus (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::minus(x,y));
	}
	const cl_UP_specialized<T> uminus (const cl_UP_specialized<T>& x)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::uminus(x));
	}
	const cl_UP_specialized<T> one ()
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::one());
	}
	const cl_UP_specialized<T> canonhom (const cl_I& x)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::canonhom(x));
	}
	const cl_UP_specialized<T> mul (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::mul(x,y));
	}
	const cl_UP_specialized<T> square (const cl_UP_specialized<T>& x)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::square(x));
	}
	const cl_UP_specialized<T> expt_pos (const cl_UP_specialized<T>& x, const cl_I& y)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::expt_pos(x,y));
	}
	const cl_UP_specialized<T> scalmul (const T& x, const cl_UP_specialized<T>& y)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::scalmul(x,y));
	}
	sintL degree (const cl_UP_specialized<T>& x)
	{
		return cl_heap_univpoly_ring::degree(x);
	}
	sintL ldegree (const cl_UP_specialized<T>& x)
	{
		return cl_heap_univpoly_ring::ldegree(x);
	}
	const cl_UP_specialized<T> monomial (const T& x, uintL e)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::monomial(cl_ring_element(cl_C_ring??,x),e));
	}
	const T coeff (const cl_UP_specialized<T>& x, uintL index)
	{
		return The(T)(cl_heap_univpoly_ring::coeff(x,index));
	}
	const cl_UP_specialized<T> create (sintL deg)
	{
		return The2(cl_UP_specialized<T>)(cl_heap_univpoly_ring::create(deg));
	}
	void set_coeff (cl_UP_specialized<T>& x, uintL index, const T& y)
	{
		cl_heap_univpoly_ring::set_coeff(x,index,cl_ring_element(cl_C_ring??,y));
	}
	void finalize (cl_UP_specialized<T>& x)
	{
		cl_heap_univpoly_ring::finalize(x);
	}
	const T eval (const cl_UP_specialized<T>& x, const T& y)
	{
		return The(T)(cl_heap_univpoly_ring::eval(x,cl_ring_element(cl_C_ring??,y)));
	}
private:
	// No need for any constructors.
	cl_heap_univpoly_specialized_ring ();
};

// Lookup of polynomial rings.
template <class T>
inline const cl_univpoly_specialized_ring<T> find_univpoly_ring (const cl_specialized_number_ring<T>& r)
{ return The(cl_univpoly_specialized_ring<T>) (find_univpoly_ring((const cl_ring&)r)); }
template <class T>
inline const cl_univpoly_specialized_ring<T> find_univpoly_ring (const cl_specialized_number_ring<T>& r, const cl_symbol& varname)
{ return The(cl_univpoly_specialized_ring<T>) (find_univpoly_ring((const cl_ring&)r,varname)); }

// Operations on polynomials.

// Add.
template <class T>
inline const cl_UP_specialized<T> operator+ (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{ return x.ring()->plus(x,y); }

// Negate.
template <class T>
inline const cl_UP_specialized<T> operator- (const cl_UP_specialized<T>& x)
	{ return x.ring()->uminus(x); }

// Subtract.
template <class T>
inline const cl_UP_specialized<T> operator- (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{ return x.ring()->minus(x,y); }

// Multiply.
template <class T>
inline const cl_UP_specialized<T> operator* (const cl_UP_specialized<T>& x, const cl_UP_specialized<T>& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
template <class T>
inline const cl_UP_specialized<T> square (const cl_UP_specialized<T>& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
template <class T>
inline const cl_UP_specialized<T> expt_pos (const cl_UP_specialized<T>& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Scalar multiplication.
// Need more discrimination on T ??
template <class T>
inline const cl_UP_specialized<T> operator* (const cl_I& x, const cl_UP_specialized<T>& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
template <class T>
inline const cl_UP_specialized<T> operator* (const cl_UP_specialized<T>& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }
template <class T>
inline const cl_UP_specialized<T> operator* (const T& x, const cl_UP_specialized<T>& y)
	{ return y.ring()->scalmul(x,y); }
template <class T>
inline const cl_UP_specialized<T> operator* (const cl_UP_specialized<T>& x, const T& y)
	{ return x.ring()->scalmul(y,x); }

// Coefficient.
template <class T>
inline const T coeff (const cl_UP_specialized<T>& x, uintL index)
	{ return x.ring()->coeff(x,index); }

// Destructive modification.
template <class T>
inline void set_coeff (cl_UP_specialized<T>& x, uintL index, const T& y)
	{ x.ring()->set_coeff(x,index,y); }
template <class T>
inline void finalize (cl_UP_specialized<T>& x)
	{ x.ring()->finalize(x); }
template <class T>
inline void cl_UP_specialized<T>::set_coeff (uintL index, const T& y)
	{ ring()->set_coeff(*this,index,y); }
template <class T>
inline void cl_UP_specialized<T>::finalize ()
	{ ring()->finalize(*this); }

// Evaluation. (No extension of the base ring allowed here for now.)
template <class T>
inline const T cl_UP_specialized<T>::operator() (const T& y) const
{
	return ring()->eval(*this,y);
}

// Derivative.
template <class T>
inline const cl_UP_specialized<T> deriv (const cl_UP_specialized<T>& x)
	{ return The(cl_UP_specialized<T>)(deriv((const cl_UP&)x)); }


#endif /* _CL_UNIVPOLY_AUX_H */
#endif

#endif /* notyet */

}  // namespace cln
