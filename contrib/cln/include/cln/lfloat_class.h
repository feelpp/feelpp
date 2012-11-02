// Concrete class of long float numbers.

#ifndef _CL_LFLOAT_CLASS_H
#define _CL_LFLOAT_CLASS_H

#include "cln/number.h"
#include "cln/float_class.h"

namespace cln {

class cl_LF : public cl_F {
public:
// Default constructor.
	cl_LF ();
// Assignment operators.
	cl_LF& operator= (const cl_LF&);
// Optimization of method pointer_p().
	bool pointer_p() const
		{ return true; }
// Faster pointer_p() gives a faster copy constructor (but not destructor!!!).
	cl_LF (const cl_LF& x);
// Other constructors.
	cl_LF (const char *);
// Private constructor.
	cl_LF (cl_private_thing);
	cl_LF (struct cl_heap_lfloat *);
// Private pointer manipulations.
	operator struct cl_heap_lfloat * () const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};
// Define this if you want the elementary cl_LF operations (+, -, *, /,
// sqrt, cl_LF_I_mul) to return results which are always the correctly
// rounded exact results, i.e. results which are correct within 0.5 ulp.
// If you don't define this, results will be correct within 0.50001 ulp,
// but often the computation will be much faster.
/* #define CL_LF_PEDANTIC */

// Private constructors.
inline cl_LF::cl_LF (cl_private_thing ptr) : cl_F (ptr) {}
// The assignment operators:
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_LF, cl_LF)
// The default constructors.
// Private pointer manipulations. Never throw away a `struct cl_heap_lfloat *'!
inline cl_LF::operator struct cl_heap_lfloat * () const
{
	struct cl_heap_lfloat * hpointer = (struct cl_heap_lfloat *) pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
extern const cl_LF cl_LF_0;
inline cl_LF::cl_LF ()
	: cl_F ((cl_private_thing) (struct cl_heap_lfloat *) cl_LF_0) {}
class cl_LF_globals_init_helper
{
	static int count;
public:
	cl_LF_globals_init_helper();
	~cl_LF_globals_init_helper();
};
static cl_LF_globals_init_helper cl_LF_globals_init_helper_instance;
#if 0 // see cl_LF_impl.h
inline cl_LF::cl_LF (struct cl_heap_lfloat * ptr)
	: cl_F ((cl_private_thing) ptr) {}
#endif
// The copy constructors.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_LF,cl_F)

}  // namespace cln

#endif /* _CL_LFLOAT_CLASS_H */
