// Ring operations.

#ifndef _CL_RING_H
#define _CL_RING_H

#include "cln/object.h"
#include "cln/malloc.h"
#include "cln/proplist.h"
#include "cln/number.h"
#include "cln/exception.h"
#include "cln/io.h"

namespace cln {

class cl_I;

// This file defines the general layout of rings, ring elements, and
// operations available on ring elements. Any subclass of `cl_ring'
// must implement these operations, with the same memory layout.
// (Because generic packages like the polynomial rings access the base
// ring's operation vectors through inline functions defined in this file.)

class cl_heap_ring;

// Rings are reference counted, but not freed immediately when they aren't
// used any more. Hence they inherit from `cl_rcpointer'.

// Vectors of function pointers are more efficient than virtual member
// functions. But it constrains us not to use multiple or virtual inheritance.
//
// Note! We are passing raw `cl_heap_ring*' pointers to the operations
// for efficiency (compared to passing `const cl_ring&', we save a memory
// access, and it is easier to cast to a `cl_heap_ring_specialized*').
// These raw pointers are meant to be used downward (in the dynamic extent
// of the call) only. If you need to save them in a data structure, cast
// to `cl_ring'; this will correctly increment the reference count.
// (This technique is safe because the inline wrapper functions make sure
// that we have a `cl_ring' somewhere containing the pointer, so there
// is no danger of dangling pointers.)
//
// Note! Because the `cl_heap_ring*' -> `cl_ring' conversion increments
// the reference count, you have to use the `cl_private_thing' -> `cl_ring'
// conversion if the reference count is already incremented.

class cl_ring : public cl_rcpointer {
public:
	// Constructor. Takes a cl_heap_ring*, increments its refcount.
	cl_ring (cl_heap_ring* r);
	// Private constructor. Doesn't increment the refcount.
	cl_ring (cl_private_thing);
	// Copy constructor.
	cl_ring (const cl_ring&);
	// Assignment operator.
	cl_ring& operator= (const cl_ring&);
	// Default constructor.
	cl_ring ();
	// Automatic dereferencing.
	cl_heap_ring* operator-> () const
	{ return (cl_heap_ring*)heappointer; }
};
CL_DEFINE_COPY_CONSTRUCTOR2(cl_ring,cl_rcpointer)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_ring,cl_ring)

// Normal constructor for `cl_ring'.
inline cl_ring::cl_ring (cl_heap_ring* r)
{ cl_inc_pointer_refcount((cl_heap*)r); pointer = r; }
// Private constructor for `cl_ring'.
inline cl_ring::cl_ring (cl_private_thing p)
{ pointer = p; }

inline bool operator== (const cl_ring& R1, const cl_ring& R2)
{ return (R1.pointer == R2.pointer); }
inline bool operator!= (const cl_ring& R1, const cl_ring& R2)
{ return (R1.pointer != R2.pointer); }
inline bool operator== (const cl_ring& R1, cl_heap_ring* R2)
{ return (R1.pointer == R2); }
inline bool operator!= (const cl_ring& R1, cl_heap_ring* R2)
{ return (R1.pointer != R2); }

// Representation of an element of a ring.
//
// In order to support true polymorphism (without C++ templates), all
// ring elements share the same basic layout:
//      cl_ring ring;     // the ring
//      cl_gcobject rep;  // representation of the element
// The representation of the element depends on the ring, of course,
// but we constrain it to be a single pointer into the heap or an immediate
// value.
//
// Any arithmetic operation on a ring R (like +, -, *) must return a value
// with ring = R. This is
// a. necessary if the computation is to proceed correctly (e.g. in cl_RA,
//    ((3/4)*4 mod 3) is 0, simplifying it to ((cl_I)4 mod (cl_I)3) = 1
//    wouldn't be correct),
// b. possible even if R is an extension ring of some ring R1 (e.g. cl_N
//    being an extension ring of cl_R). Automatic retraction from R to R1
//    can be done through dynamic typing: An element of R which happens
//    to lie in R1 is stored using the internal representation of R1,
//    but with ring = R. Elements of R1 and R\R1 can be distinguished
//    through rep's type.
// c. an advantage for the implementation of polynomials and other
//    entities which contain many elements of the same ring. They need
//    to store only the elements' representations, and a single pointer
//    to the ring.
//
// The ring operations exist in two versions:
// - Low-level version, which only operates on the representation.
// - High-level version, which operates on full cl_ring_elements.
// We make this distinction for performance: Multiplication of polynomials
// over Z/nZ, operating on the high-level operations, spends 40% of its
// computing time with packing and unpacking of cl_ring_elements.
// The low-level versions have an underscore prepended and are unsafe.

class _cl_ring_element {
public:
	cl_gcobject rep;	// representation of the element
	// Default constructor.
	_cl_ring_element ();
public: /* ugh */
	// Constructor.
	_cl_ring_element (const cl_heap_ring* R, const cl_gcobject& r) : rep (as_cl_private_thing(r)) { (void)R; }
	_cl_ring_element (const cl_ring& R, const cl_gcobject& r) : rep (as_cl_private_thing(r)) { (void)R; }
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

class cl_ring_element : public _cl_ring_element {
protected:
	cl_ring _ring;			// ring
public:
	const cl_ring& ring () const { return _ring; }
	// Default constructor.
	cl_ring_element ();
public: /* ugh */
	// Constructor.
	cl_ring_element (const cl_ring& R, const cl_gcobject& r) : _cl_ring_element (R,r), _ring (R) {}
	cl_ring_element (const cl_ring& R, const _cl_ring_element& r) : _cl_ring_element (r), _ring (R) {}
public:	// Debugging output.
	void debug_print () const;
	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

// The ring operations are encoded as vectors of function pointers. You
// can add more operations to the end of each vector or add new vectors,
// but you must not reorder the operations nor reorder the vectors nor
// change the functions' signatures incompatibly.

// There should ideally be a template class for each vector, but unfortunately
// you lose the ability to initialize the vector using "= { ... }" syntax
// when you subclass it.

struct _cl_ring_setops {
	// print
	void (* fprint) (cl_heap_ring* R, std::ostream& stream, const _cl_ring_element& x);
	// equality
	bool (* equal) (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y);
	// ...
};
struct _cl_ring_addops {
	// 0
	const _cl_ring_element (* zero) (cl_heap_ring* R);
	bool (* zerop) (cl_heap_ring* R, const _cl_ring_element& x);
	// x+y
	const _cl_ring_element (* plus) (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y);
	// x-y
	const _cl_ring_element (* minus) (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y);
	// -x
	const _cl_ring_element (* uminus) (cl_heap_ring* R, const _cl_ring_element& x);
	// ...
};
struct _cl_ring_mulops {
	// 1
	const _cl_ring_element (* one) (cl_heap_ring* R);
	// canonical homomorphism
	const _cl_ring_element (* canonhom) (cl_heap_ring* R, const cl_I& x);
	// x*y
	const _cl_ring_element (* mul) (cl_heap_ring* R, const _cl_ring_element& x, const _cl_ring_element& y);
	// x^2
	const _cl_ring_element (* square) (cl_heap_ring* R, const _cl_ring_element& x);
	// x^y, y Integer >0
	const _cl_ring_element (* expt_pos) (cl_heap_ring* R, const _cl_ring_element& x, const cl_I& y);
	// ...
};
  typedef const _cl_ring_setops  cl_ring_setops;
  typedef const _cl_ring_addops  cl_ring_addops;
  typedef const _cl_ring_mulops  cl_ring_mulops;

// Representation of a ring in memory.

class cl_heap_ring : public cl_heap {
public:
	// Allocation.
	void* operator new (size_t size) { return malloc_hook(size); }
	// Deallocation.
	void operator delete (void* ptr) { free_hook(ptr); }
private:
	cl_property_list properties;
protected:
	cl_ring_setops* setops;
	cl_ring_addops* addops;
	cl_ring_mulops* mulops;
public:
	// More information comes here.
	// ...
public:
	// Low-level operations.
	void _fprint (std::ostream& stream, const _cl_ring_element& x)
		{ setops->fprint(this,stream,x); }
	bool _equal (const _cl_ring_element& x, const _cl_ring_element& y)
		{ return setops->equal(this,x,y); }
	const _cl_ring_element _zero ()
		{ return addops->zero(this); }
	bool _zerop (const _cl_ring_element& x)
		{ return addops->zerop(this,x); }
	const _cl_ring_element _plus (const _cl_ring_element& x, const _cl_ring_element& y)
		{ return addops->plus(this,x,y); }
	const _cl_ring_element _minus (const _cl_ring_element& x, const _cl_ring_element& y)
		{ return addops->minus(this,x,y); }
	const _cl_ring_element _uminus (const _cl_ring_element& x)
		{ return addops->uminus(this,x); }
	const _cl_ring_element _one ()
		{ return mulops->one(this); }
	const _cl_ring_element _canonhom (const cl_I& x)
		{ return mulops->canonhom(this,x); }
	const _cl_ring_element _mul (const _cl_ring_element& x, const _cl_ring_element& y)
		{ return mulops->mul(this,x,y); }
	const _cl_ring_element _square (const _cl_ring_element& x)
		{ return mulops->square(this,x); }
	const _cl_ring_element _expt_pos (const _cl_ring_element& x, const cl_I& y)
		{ return mulops->expt_pos(this,x,y); }
	// High-level operations.
	void fprint (std::ostream& stream, const cl_ring_element& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		_fprint(stream,x);
	}
	bool equal (const cl_ring_element& x, const cl_ring_element& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return _equal(x,y);
	}
	const cl_ring_element zero ()
	{
		return cl_ring_element(this,_zero());
	}
	bool zerop (const cl_ring_element& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _zerop(x);
	}
	const cl_ring_element plus (const cl_ring_element& x, const cl_ring_element& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_ring_element(this,_plus(x,y));
	}
	const cl_ring_element minus (const cl_ring_element& x, const cl_ring_element& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_ring_element(this,_minus(x,y));
	}
	const cl_ring_element uminus (const cl_ring_element& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_ring_element(this,_uminus(x));
	}
	const cl_ring_element one ()
	{
		return cl_ring_element(this,_one());
	}
	const cl_ring_element canonhom (const cl_I& x)
	{
		return cl_ring_element(this,_canonhom(x));
	}
	const cl_ring_element mul (const cl_ring_element& x, const cl_ring_element& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_ring_element(this,_mul(x,y));
	}
	const cl_ring_element square (const cl_ring_element& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_ring_element(this,_square(x));
	}
	const cl_ring_element expt_pos (const cl_ring_element& x, const cl_I& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_ring_element(this,_expt_pos(x,y));
	}
	// Property operations.
	cl_property* get_property (const cl_symbol& key)
		{ return properties.get_property(key); }
	void add_property (cl_property* new_property)
		{ properties.add_property(new_property); }
// Constructor.
	cl_heap_ring (cl_ring_setops* setopv, cl_ring_addops* addopv, cl_ring_mulops* mulopv)
		: setops (setopv), addops (addopv), mulops (mulopv)
		{ refcount = 0; } // will be incremented by the `cl_ring' constructor
};
#define SUBCLASS_cl_heap_ring() \
public:									  \
	/* Allocation. */						  \
	void* operator new (size_t size) { return malloc_hook(size); } \
	/* Deallocation. */						  \
	void operator delete (void* ptr) { free_hook(ptr); }

// Operations on ring elements.

// Output.
inline void fprint (std::ostream& stream, const cl_ring_element& x)
	{ x.ring()->fprint(stream,x); }
CL_DEFINE_PRINT_OPERATOR(cl_ring_element)

// Add.
inline const cl_ring_element operator+ (const cl_ring_element& x, const cl_ring_element& y)
	{ return x.ring()->plus(x,y); }

// Negate.
inline const cl_ring_element operator- (const cl_ring_element& x)
	{ return x.ring()->uminus(x); }

// Subtract.
inline const cl_ring_element operator- (const cl_ring_element& x, const cl_ring_element& y)
	{ return x.ring()->minus(x,y); }

// Equality.
inline bool operator== (const cl_ring_element& x, const cl_ring_element& y)
	{ return x.ring()->equal(x,y); }
inline bool operator!= (const cl_ring_element& x, const cl_ring_element& y)
	{ return !x.ring()->equal(x,y); }

// Compare against 0.
inline bool zerop (const cl_ring_element& x)
	{ return x.ring()->zerop(x); }

// Multiply.
inline const cl_ring_element operator* (const cl_ring_element& x, const cl_ring_element& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
inline const cl_ring_element square (const cl_ring_element& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
inline const cl_ring_element expt_pos (const cl_ring_element& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Scalar multiplication.
// [Is this operation worth being specially optimized for the case of
// polynomials?? Polynomials have a faster scalar multiplication.
// We should use it.??]
inline const cl_ring_element operator* (const cl_I& x, const cl_ring_element& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
inline const cl_ring_element operator* (const cl_ring_element& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }


// Ring of uninitialized elements.
// Any operation results in an exception being thrown.

// Thrown when an attempt is made to perform an operation on an uninitialized ring.
class uninitialized_ring_exception : public runtime_exception {
public:
	uninitialized_ring_exception ();
};

// Thrown when a ring element is uninitialized.
class uninitialized_exception : public runtime_exception {
public:
	explicit uninitialized_exception (const _cl_ring_element& obj);
	uninitialized_exception (const _cl_ring_element& obj_x, const _cl_ring_element& obj_y);
};

extern const cl_ring cl_no_ring;
extern cl_class cl_class_no_ring;

class cl_no_ring_init_helper
{
	static int count;
public:
	cl_no_ring_init_helper();
	~cl_no_ring_init_helper();
};
static cl_no_ring_init_helper cl_no_ring_init_helper_instance;

inline cl_ring::cl_ring ()
	: cl_rcpointer (as_cl_private_thing(cl_no_ring)) {}
inline _cl_ring_element::_cl_ring_element ()
	: rep ((cl_private_thing) cl_combine(cl_FN_tag,0)) {}
inline cl_ring_element::cl_ring_element ()
	: _cl_ring_element (), _ring () {}


// Support for built-in number rings.
// Beware, they are not optimally efficient.

template <class T>
struct cl_number_ring_ops {
	bool (* contains) (const cl_number&);
	bool (* equal) (const T&, const T&);
	bool (* zerop) (const T&);
	const T (* plus) (const T&, const T&);
	const T (* minus) (const T&, const T&);
	const T (* uminus) (const T&);
	const T (* mul) (const T&, const T&);
	const T (* square) (const T&);
	const T (* expt_pos) (const T&, const cl_I&);
};
class cl_heap_number_ring : public cl_heap_ring {
public:
	cl_number_ring_ops<cl_number>* ops;
	// Constructor.
	cl_heap_number_ring (cl_ring_setops* setopv, cl_ring_addops* addopv, cl_ring_mulops* mulopv, cl_number_ring_ops<cl_number>* opv)
		: cl_heap_ring (setopv,addopv,mulopv), ops (opv) {}
};

class cl_number_ring : public cl_ring {
public:
	cl_number_ring (cl_heap_number_ring* r)
		: cl_ring (r) {}
};

template <class T>
class cl_specialized_number_ring : public cl_number_ring {
public:
	cl_specialized_number_ring ();
};

// Type test.
inline bool instanceof (const cl_number& x, const cl_number_ring& R)
{
	return ((cl_heap_number_ring*) R.heappointer)->ops->contains(x);
}


// Hack section.

// Conversions to subtypes without checking:
// The2(cl_MI)(x) converts x to a cl_MI, without change of representation!
  #define The(type)  *(const type *) & cl_identity
  #define The2(type)  *(const type *) & cl_identity2
// This inline function is for type checking purposes only.
  inline const cl_ring& cl_identity (const cl_ring& r) { return r; }
  inline const cl_ring_element& cl_identity2 (const cl_ring_element& x) { return x; }
  inline const cl_gcobject& cl_identity (const _cl_ring_element& x) { return x.rep; }


// Debugging support.
#ifdef CL_DEBUG
extern int cl_ring_debug_module;
CL_FORCE_LINK(cl_ring_debug_dummy, cl_ring_debug_module)
#endif

}  // namespace cln

#endif /* _CL_RING_H */
