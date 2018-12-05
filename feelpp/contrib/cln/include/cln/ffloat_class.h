// Concrete class of single float numbers.

#ifndef _CL_FFLOAT_CLASS_H
#define _CL_FFLOAT_CLASS_H

#include "cln/number.h"
#include "cln/float_class.h"

namespace cln {

class cl_FF : public cl_F {
public:
// Default constructor.
	cl_FF ();
// Assignment operators.
	cl_FF& operator= (const cl_FF&);
// Optimization of method pointer_p().
	bool pointer_p() const
#if defined(CL_WIDE_POINTERS)
		{ return false; }
#else
		{ return true; }
#endif
// Faster pointer_p() gives a faster copy constructor (but not destructor!!!).
	cl_FF (const cl_FF& x);
// Constructors and assignment operators from C numeric types.
	cl_FF (const float);
	cl_FF& operator= (const float);
// Other constructors.
	cl_FF (const char *);
// Private constructor.
	cl_FF (cl_private_thing);
#if defined(CL_WIDE_POINTERS)
	cl_FF (struct cl_heap_ffloat * /* NULL! */, cl_uint);
#else
	cl_FF (struct cl_heap_ffloat *);
// Private pointer manipulations.
	operator struct cl_heap_ffloat * () const;
#endif
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

// Private constructors.
inline cl_FF::cl_FF (cl_private_thing ptr) : cl_F (ptr) {}
// The assignment operators:
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_FF, cl_FF)
// The default constructors.
#if defined(CL_WIDE_POINTERS)
inline cl_FF::cl_FF ()
	: cl_F ((cl_private_thing) cl_combine(cl_FF_tag,0)) {}
#else
// Private pointer manipulations. Never throw away a `struct cl_heap_ffloat *'!
inline cl_FF::operator struct cl_heap_ffloat * () const
{
	struct cl_heap_ffloat * hpointer = (struct cl_heap_ffloat *) pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
extern const cl_FF cl_FF_0;
inline cl_FF::cl_FF ()
	: cl_F ((cl_private_thing) (struct cl_heap_ffloat *) cl_FF_0) {}
class cl_FF_globals_init_helper
{
	static int count;
public:
	cl_FF_globals_init_helper();
	~cl_FF_globals_init_helper();
};
static cl_FF_globals_init_helper cl_FF_globals_init_helper_instance;
#if 0 // see cl_FF.h
inline cl_FF::cl_FF (struct cl_heap_ffloat * ptr)
	: cl_F ((cl_private_thing) ptr) {}
#endif
#endif
// The copy constructors.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_FF,cl_F)
// Constructors and assignment operators from C numeric types.
CL_DEFINE_FLOAT_CONSTRUCTOR(cl_FF)

}  // namespace cln

#endif /* _CL_FFLOAT_CLASS_H */
