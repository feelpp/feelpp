// Concrete class of short float numbers.

#ifndef _CL_SFLOAT_CLASS_H
#define _CL_SFLOAT_CLASS_H

#include "cln/number.h"
#include "cln/float_class.h"

namespace cln {

class cl_SF : public cl_F {
public:
// Default constructor.
	cl_SF ();
// Assignment operators.
	cl_SF& operator= (const cl_SF&);
// Optimization of method pointer_p().
	bool pointer_p() const
		{ return false; }
// Faster pointer_p() gives a faster copy constructor (but not destructor!!!).
	cl_SF (const cl_SF& x);
// Other constructors.
	cl_SF (const char *);
// Private constructor.
	cl_SF (cl_private_thing);
	cl_SF (struct cl_sfloat * /* NULL! */, cl_uint);
public:	// Ability to place an object at a given address.
	void* operator new (size_t size) { return malloc_hook(size); }
	void* operator new (size_t size, void* ptr) { (void)size; return ptr; }
	void operator delete (void* ptr) { free_hook(ptr); }
};

// Private constructors.
inline cl_SF::cl_SF (cl_private_thing ptr) : cl_F (ptr) {}
// The assignment operators:
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SF, cl_SF)
// The default constructors.
inline cl_SF::cl_SF ()
	: cl_F ((cl_private_thing) cl_combine(cl_SF_tag,0)) {}
// The copy constructors.
CL_DEFINE_COPY_CONSTRUCTOR2(cl_SF,cl_F)

}  // namespace cln

#endif /* _CL_SFLOAT_CLASS_H */
