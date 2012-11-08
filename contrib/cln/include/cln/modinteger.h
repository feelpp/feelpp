// Modular integer operations.

#ifndef _CL_MODINTEGER_H
#define _CL_MODINTEGER_H

#include "cln/object.h"
#include "cln/ring.h"
#include "cln/integer.h"
#include "cln/random.h"
#include "cln/malloc.h"
#include "cln/io.h"
#include "cln/proplist.h"
#include "cln/condition.h"
#include "cln/exception.h"
#undef random // Linux defines random() as a macro!

namespace cln {

// Representation of an element of a ring Z/mZ.

// To protect against mixing elements of different modular rings, such as
// (3 mod 4) + (2 mod 5), every modular integer carries its ring in itself.


// Representation of a ring Z/mZ.

class cl_heap_modint_ring;

class cl_modint_ring : public cl_ring {
public:
	// Default constructor.
	cl_modint_ring ();
	// Constructor. Takes a cl_heap_modint_ring*, increments its refcount.
	cl_modint_ring (cl_heap_modint_ring* r);
	// Copy constructor.
	cl_modint_ring (const cl_modint_ring&);
	// Assignment operator.
	cl_modint_ring& operator= (const cl_modint_ring&);
	// Automatic dereferencing.
	cl_heap_modint_ring* operator-> () const
	{ return (cl_heap_modint_ring*)heappointer; }
};

// Z/0Z
extern const cl_modint_ring cl_modint0_ring;
// Default constructor. This avoids dealing with NULL pointers.
inline cl_modint_ring::cl_modint_ring ()
	: cl_ring (as_cl_private_thing(cl_modint0_ring)) {}

class cl_MI_init_helper
{
	static int count;
public:
	cl_MI_init_helper();
	~cl_MI_init_helper();
};
static cl_MI_init_helper cl_MI_init_helper_instance;

// Copy constructor and assignment operator.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_modint_ring,cl_ring)
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_modint_ring,cl_modint_ring)

// Normal constructor for `cl_modint_ring'.
inline cl_modint_ring::cl_modint_ring (cl_heap_modint_ring* r)
	: cl_ring ((cl_private_thing) (cl_inc_pointer_refcount((cl_heap*)r), r)) {}

// Operations on modular integer rings.

inline bool operator== (const cl_modint_ring& R1, const cl_modint_ring& R2)
{ return (R1.pointer == R2.pointer); }
inline bool operator!= (const cl_modint_ring& R1, const cl_modint_ring& R2)
{ return (R1.pointer != R2.pointer); }
inline bool operator== (const cl_modint_ring& R1, cl_heap_modint_ring* R2)
{ return (R1.pointer == R2); }
inline bool operator!= (const cl_modint_ring& R1, cl_heap_modint_ring* R2)
{ return (R1.pointer != R2); }


// Condition raised when a probable prime is discovered to be composite.
struct cl_composite_condition : public cl_condition {
	SUBCLASS_cl_condition()
	cl_I p; // the non-prime
	cl_I factor; // a nontrivial factor, or 0
	// Constructors.
	cl_composite_condition (const cl_I& _p)
		: p (_p), factor (0)
		{ print(std::cerr); }
	cl_composite_condition (const cl_I& _p, const cl_I& _f)
		: p (_p), factor (_f)
		{ print(std::cerr); }
	// Implement general condition methods.
	const char * name () const;
	void print (std::ostream&) const;
	~cl_composite_condition () {}
};


// Representation of an element of a ring Z/mZ.

class _cl_MI /* cf. _cl_ring_element */ {
public:
	cl_I rep;		// representative, integer >=0, <m
				// (maybe the Montgomery representative!)
	// Default constructor.
	_cl_MI () : rep () {}
public: /* ugh */
	// Constructor.
	_cl_MI (const cl_heap_modint_ring* R, const cl_I& r) : rep (r) { (void)R; }
	_cl_MI (const cl_modint_ring& R, const cl_I& r) : rep (r) { (void)R; }
public:
	// Conversion.
	CL_DEFINE_CONVERTER(_cl_ring_element)
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

class cl_MI /* cf. cl_ring_element */ : public _cl_MI {
protected:
	cl_modint_ring _ring;	// ring Z/mZ
public:
	const cl_modint_ring& ring () const { return _ring; }
	// Default constructor.
	cl_MI () : _cl_MI (), _ring () {}
public: /* ugh */
	// Constructor.
	cl_MI (const cl_modint_ring& R, const cl_I& r) : _cl_MI (R,r), _ring (R) {}
	cl_MI (const cl_modint_ring& R, const _cl_MI& r) : _cl_MI (r), _ring (R) {}
public:
	// Conversion.
	CL_DEFINE_CONVERTER(cl_ring_element)
	// Debugging output.
	void debug_print () const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};


// Representation of an element of a ring Z/mZ or an exception.

class cl_MI_x {
private:
	cl_MI value;
public:
	cl_composite_condition* condition;
	// Constructors.
	cl_MI_x (cl_composite_condition* c) : value (), condition (c) {}
	cl_MI_x (const cl_MI& x) : value (x), condition (NULL) {}
	// Cast operators.
	//operator cl_MI& () { if (condition) throw runtime_exception(); return value; }
	//operator const cl_MI& () const { if (condition) throw runtime_exception(); return value; }
	operator cl_MI () const { if (condition) throw runtime_exception(); return value; }
};


// Ring operations.

struct _cl_modint_setops /* cf. _cl_ring_setops */ {
	// print
	void (* fprint) (cl_heap_modint_ring* R, std::ostream& stream, const _cl_MI& x);
	// equality
	bool (* equal) (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y);
	// random number
	const _cl_MI (* random) (cl_heap_modint_ring* R, random_state& randomstate);
};
struct _cl_modint_addops /* cf. _cl_ring_addops */ {
	// 0
	const _cl_MI (* zero) (cl_heap_modint_ring* R);
	bool (* zerop) (cl_heap_modint_ring* R, const _cl_MI& x);
	// x+y
	const _cl_MI (* plus) (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y);
	// x-y
	const _cl_MI (* minus) (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y);
	// -x
	const _cl_MI (* uminus) (cl_heap_modint_ring* R, const _cl_MI& x);
};
struct _cl_modint_mulops /* cf. _cl_ring_mulops */ {
	// 1
	const _cl_MI (* one) (cl_heap_modint_ring* R);
	// canonical homomorphism
	const _cl_MI (* canonhom) (cl_heap_modint_ring* R, const cl_I& x);
	// x*y
	const _cl_MI (* mul) (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y);
	// x^2
	const _cl_MI (* square) (cl_heap_modint_ring* R, const _cl_MI& x);
	// x^y, y Integer >0
	const _cl_MI (* expt_pos) (cl_heap_modint_ring* R, const _cl_MI& x, const cl_I& y);
	// x^-1
	const cl_MI_x (* recip) (cl_heap_modint_ring* R, const _cl_MI& x);
	// x*y^-1
	const cl_MI_x (* div) (cl_heap_modint_ring* R, const _cl_MI& x, const _cl_MI& y);
	// x^y, y Integer
	const cl_MI_x (* expt) (cl_heap_modint_ring* R, const _cl_MI& x, const cl_I& y);
	// x -> x mod m for x>=0
	const cl_I (* reduce_modulo) (cl_heap_modint_ring* R, const cl_I& x);
	// some inverse of canonical homomorphism
	const cl_I (* retract) (cl_heap_modint_ring* R, const _cl_MI& x);
};
  typedef const _cl_modint_setops  cl_modint_setops;
  typedef const _cl_modint_addops  cl_modint_addops;
  typedef const _cl_modint_mulops  cl_modint_mulops;

// Representation of the ring Z/mZ.

// Currently rings are garbage collected only when they are not referenced
// any more and when the ring table gets full.

// Modular integer rings are kept unique in memory. This way, ring equality
// can be checked very efficiently by a simple pointer comparison.

class cl_heap_modint_ring /* cf. cl_heap_ring */ : public cl_heap {
	SUBCLASS_cl_heap_ring()
private:
	cl_property_list properties;
protected:
	cl_modint_setops* setops;
	cl_modint_addops* addops;
	cl_modint_mulops* mulops;
public:
	cl_I modulus;	// m, normalized to be >= 0
public:
	// Low-level operations.
	void _fprint (std::ostream& stream, const _cl_MI& x)
		{ setops->fprint(this,stream,x); }
	bool _equal (const _cl_MI& x, const _cl_MI& y)
		{ return setops->equal(this,x,y); }
	const _cl_MI _random (random_state& randomstate)
		{ return setops->random(this,randomstate); }
	const _cl_MI _zero ()
		{ return addops->zero(this); }
	bool _zerop (const _cl_MI& x)
		{ return addops->zerop(this,x); }
	const _cl_MI _plus (const _cl_MI& x, const _cl_MI& y)
		{ return addops->plus(this,x,y); }
	const _cl_MI _minus (const _cl_MI& x, const _cl_MI& y)
		{ return addops->minus(this,x,y); }
	const _cl_MI _uminus (const _cl_MI& x)
		{ return addops->uminus(this,x); }
	const _cl_MI _one ()
		{ return mulops->one(this); }
	const _cl_MI _canonhom (const cl_I& x)
		{ return mulops->canonhom(this,x); }
	const _cl_MI _mul (const _cl_MI& x, const _cl_MI& y)
		{ return mulops->mul(this,x,y); }
	const _cl_MI _square (const _cl_MI& x)
		{ return mulops->square(this,x); }
	const _cl_MI _expt_pos (const _cl_MI& x, const cl_I& y)
		{ return mulops->expt_pos(this,x,y); }
	const cl_MI_x _recip (const _cl_MI& x)
		{ return mulops->recip(this,x); }
	const cl_MI_x _div (const _cl_MI& x, const _cl_MI& y)
		{ return mulops->div(this,x,y); }
	const cl_MI_x _expt (const _cl_MI& x, const cl_I& y)
		{ return mulops->expt(this,x,y); }
	const cl_I _reduce_modulo (const cl_I& x)
		{ return mulops->reduce_modulo(this,x); }
	const cl_I _retract (const _cl_MI& x)
		{ return mulops->retract(this,x); }
	// High-level operations.
	void fprint (std::ostream& stream, const cl_MI& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		_fprint(stream,x);
	}
	bool equal (const cl_MI& x, const cl_MI& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return _equal(x,y);
	}
	const cl_MI random (random_state& randomstate = default_random_state)
	{
		return cl_MI(this,_random(randomstate));
	}
	const cl_MI zero ()
	{
		return cl_MI(this,_zero());
	}
	bool zerop (const cl_MI& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _zerop(x);
	}
	const cl_MI plus (const cl_MI& x, const cl_MI& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_MI(this,_plus(x,y));
	}
	const cl_MI minus (const cl_MI& x, const cl_MI& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_MI(this,_minus(x,y));
	}
	const cl_MI uminus (const cl_MI& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_MI(this,_uminus(x));
	}
	const cl_MI one ()
	{
		return cl_MI(this,_one());
	}
	const cl_MI canonhom (const cl_I& x)
	{
		return cl_MI(this,_canonhom(x));
	}
	const cl_MI mul (const cl_MI& x, const cl_MI& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return cl_MI(this,_mul(x,y));
	}
	const cl_MI square (const cl_MI& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_MI(this,_square(x));
	}
	const cl_MI expt_pos (const cl_MI& x, const cl_I& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return cl_MI(this,_expt_pos(x,y));
	}
	const cl_MI_x recip (const cl_MI& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _recip(x);
	}
	const cl_MI_x div (const cl_MI& x, const cl_MI& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		if (!(y.ring() == this)) throw runtime_exception();
		return _div(x,y);
	}
	const cl_MI_x expt (const cl_MI& x, const cl_I& y)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _expt(x,y);
	}
	const cl_I reduce_modulo (const cl_I& x)
	{
		return _reduce_modulo(x);
	}
	const cl_I retract (const cl_MI& x)
	{
		if (!(x.ring() == this)) throw runtime_exception();
		return _retract(x);
	}
	// Miscellaneous.
	sintC bits; // number of bits needed to represent a representative, or -1
	int log2_bits; // log_2(bits), or -1
	// Property operations.
	cl_property* get_property (const cl_symbol& key)
		{ return properties.get_property(key); }
	void add_property (cl_property* new_property)
		{ properties.add_property(new_property); }
// Constructor / destructor.
	cl_heap_modint_ring (cl_I m, cl_modint_setops*, cl_modint_addops*, cl_modint_mulops*);
	~cl_heap_modint_ring () {}
};
#define SUBCLASS_cl_heap_modint_ring() \
  SUBCLASS_cl_heap_ring()

// Lookup or create a modular integer ring  Z/mZ
extern const cl_modint_ring find_modint_ring (const cl_I& m);
static cl_MI_init_helper cl_MI_init_helper_instance2;

// Operations on modular integers.

// Output.
inline void fprint (std::ostream& stream, const cl_MI& x)
	{ x.ring()->fprint(stream,x); }
CL_DEFINE_PRINT_OPERATOR(cl_MI)

// Add.
inline const cl_MI operator+ (const cl_MI& x, const cl_MI& y)
	{ return x.ring()->plus(x,y); }
inline const cl_MI operator+ (const cl_MI& x, const cl_I& y)
	{ return x.ring()->plus(x,x.ring()->canonhom(y)); }
inline const cl_MI operator+ (const cl_I& x, const cl_MI& y)
	{ return y.ring()->plus(y.ring()->canonhom(x),y); }

// Negate.
inline const cl_MI operator- (const cl_MI& x)
	{ return x.ring()->uminus(x); }

// Subtract.
inline const cl_MI operator- (const cl_MI& x, const cl_MI& y)
	{ return x.ring()->minus(x,y); }
inline const cl_MI operator- (const cl_MI& x, const cl_I& y)
	{ return x.ring()->minus(x,x.ring()->canonhom(y)); }
inline const cl_MI operator- (const cl_I& x, const cl_MI& y)
	{ return y.ring()->minus(y.ring()->canonhom(x),y); }

// Shifts.
extern const cl_MI operator<< (const cl_MI& x, sintC y); // assume 0 <= y < 2^(intCsize-1)
extern const cl_MI operator>> (const cl_MI& x, sintC y); // assume m odd, 0 <= y < 2^(intCsize-1)

// Equality.
inline bool operator== (const cl_MI& x, const cl_MI& y)
	{ return x.ring()->equal(x,y); }
inline bool operator!= (const cl_MI& x, const cl_MI& y)
	{ return !x.ring()->equal(x,y); }
inline bool operator== (const cl_MI& x, const cl_I& y)
	{ return x.ring()->equal(x,x.ring()->canonhom(y)); }
inline bool operator!= (const cl_MI& x, const cl_I& y)
	{ return !x.ring()->equal(x,x.ring()->canonhom(y)); }
inline bool operator== (const cl_I& x, const cl_MI& y)
	{ return y.ring()->equal(y.ring()->canonhom(x),y); }
inline bool operator!= (const cl_I& x, const cl_MI& y)
	{ return !y.ring()->equal(y.ring()->canonhom(x),y); }

// Compare against 0.
inline bool zerop (const cl_MI& x)
	{ return x.ring()->zerop(x); }

// Multiply.
inline const cl_MI operator* (const cl_MI& x, const cl_MI& y)
	{ return x.ring()->mul(x,y); }

// Squaring.
inline const cl_MI square (const cl_MI& x)
	{ return x.ring()->square(x); }

// Exponentiation x^y, where y > 0.
inline const cl_MI expt_pos (const cl_MI& x, const cl_I& y)
	{ return x.ring()->expt_pos(x,y); }

// Reciprocal.
inline const cl_MI recip (const cl_MI& x)
	{ return x.ring()->recip(x); }

// Division.
inline const cl_MI div (const cl_MI& x, const cl_MI& y)
	{ return x.ring()->div(x,y); }
inline const cl_MI div (const cl_MI& x, const cl_I& y)
	{ return x.ring()->div(x,x.ring()->canonhom(y)); }
inline const cl_MI div (const cl_I& x, const cl_MI& y)
	{ return y.ring()->div(y.ring()->canonhom(x),y); }

// Exponentiation x^y.
inline const cl_MI expt (const cl_MI& x, const cl_I& y)
	{ return x.ring()->expt(x,y); }

// Scalar multiplication.
inline const cl_MI operator* (const cl_I& x, const cl_MI& y)
	{ return y.ring()->mul(y.ring()->canonhom(x),y); }
inline const cl_MI operator* (const cl_MI& x, const cl_I& y)
	{ return x.ring()->mul(x.ring()->canonhom(y),x); }

// TODO: implement gcd, index (= gcd), unitp, sqrtp


// Debugging support.
#ifdef CL_DEBUG
extern int cl_MI_debug_module;
CL_FORCE_LINK(cl_MI_debug_dummy, cl_MI_debug_module)
#endif

}  // namespace cln

#endif /* _CL_MODINTEGER_H */
