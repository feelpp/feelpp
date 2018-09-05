// Concrete class of double float numbers.

#ifndef _CL_DFLOAT_CLASS_H
#define _CL_DFLOAT_CLASS_H

#include "cln/number.h"
#include "cln/float_class.h"

namespace cln {

class cl_DF : public cl_F {
public:
// Default constructor.
	cl_DF ();
// Assignment operators.
	cl_DF& operator= (const cl_DF&);
// Optimization of method pointer_p().
	bool pointer_p() const
		{ return true; }
// Faster pointer_p() gives a faster copy constructor (but not destructor!!!).
	cl_DF (const cl_DF& x);
// Constructors and assignment operators from C numeric types.
	cl_DF (const double);
	cl_DF& operator= (const double);
// Other constructors.
	cl_DF (const char *);
// Private constructor.
	cl_DF (cl_private_thing);
	cl_DF (struct cl_heap_dfloat *);
// Private pointer manipulations.
	operator struct cl_heap_dfloat * () const;
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
private:
// Friend declarations. They are for the compiler. Just ignore them.
};

// Private constructors.
inline cl_DF::cl_DF (cl_private_thing ptr) : cl_F (ptr) {}
// The assignment operators:
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_DF, cl_DF)
// The default constructors.
// Private pointer manipulations. Never throw away a `struct cl_heap_dfloat *'!
inline cl_DF::operator struct cl_heap_dfloat * () const
{
	struct cl_heap_dfloat * hpointer = (struct cl_heap_dfloat *) pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
extern const cl_DF cl_DF_0;
inline cl_DF::cl_DF ()
	: cl_F ((cl_private_thing) (struct cl_heap_dfloat *) cl_DF_0) {}
class cl_DF_globals_init_helper
{
	static int count;
public:
	cl_DF_globals_init_helper();
	~cl_DF_globals_init_helper();
};
static cl_DF_globals_init_helper cl_DF_globals_init_helper_instance;

#if 0 // see cl_DF.h
inline cl_DF::cl_DF (struct cl_heap_dfloat * ptr)
	: cl_F ((cl_private_thing) ptr) {}
#endif
// The copy constructors.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_DF,cl_F)
// Constructors and assignment operators from C numeric types.
CL_DEFINE_DOUBLE_CONSTRUCTOR(cl_DF)

}  // namespace cln

#endif /* _CL_DFLOAT_CLASS_H */
