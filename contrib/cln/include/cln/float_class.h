// Abstract class of floating-point numbers.

#ifndef _CL_FLOAT_CLASS_H
#define _CL_FLOAT_CLASS_H

#include "cln/number.h"
#include "cln/real_class.h"

namespace cln {

class cl_F : public cl_R {
public:
// Default constructor.
	cl_F ();
// Copy constructor.
	cl_F (const cl_F&);
// Converters.
// Assignment operators.
	cl_F& operator= (const cl_F&);
// Constructors and assignment operators from C numeric types.
	cl_F (const float);
	cl_F (const double);
	cl_F& operator= (const float);
	cl_F& operator= (const double);
// Other constructors.
	cl_F (const char *);
// Private constructor.
	cl_F (cl_private_thing);
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
private:
// Friend declarations. They are for the compiler. Just ignore them.
};

// Private constructors.
inline cl_F::cl_F (cl_private_thing ptr) : cl_R (ptr) {}
// The assignment operators:
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_F, cl_F)
// The default constructors.
inline cl_F::cl_F ()
	: cl_R ((cl_private_thing) cl_combine(cl_SF_tag,0)) {}
// The copy constructors.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_F,cl_R)
// Constructors and assignment operators from C numeric types.
CL_DEFINE_FLOAT_CONSTRUCTOR(cl_F)
CL_DEFINE_DOUBLE_CONSTRUCTOR(cl_F)

}  // namespace cln

#endif /* _CL_FLOAT_CLASS_H */
