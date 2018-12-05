// Simple vectors of numbers.

#ifndef _CL_SV_NUMBER_H
#define _CL_SV_NUMBER_H

#include "cln/number.h"
#include "cln/SV.h"
#include "cln/io.h"

namespace cln {

typedef cl_heap_SV<cl_number> cl_heap_SV_number;

struct cl_SV_number : public cl_SV<cl_number,cl_SV_any> {
public:
	// Constructors.
	cl_SV_number ();
	cl_SV_number (const cl_SV_number&);
	explicit cl_SV_number (std::size_t len);
	// Assignment operators.
	cl_SV_number& operator= (const cl_SV_number&);
	// Private pointer manipulations.
	operator cl_heap_SV_number* () const;
	cl_SV_number (cl_heap_SV_number* p) : cl_SV<cl_number,cl_SV_any> (p) {}
	cl_SV_number (cl_private_thing p) : cl_SV<cl_number,cl_SV_any> (p) {}
};
inline cl_SV_number::cl_SV_number (const cl_SV_number& x) : cl_SV<cl_number,cl_SV_any> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SV_number,cl_SV_number)
// Returns a new simple vector with uninitialized contents.
extern cl_heap_SV_number* cl_make_heap_SV_number_uninit (std::size_t len);
// Returns a new simple vector with all elements initialized to 0.
extern cl_heap_SV_number* cl_make_heap_SV_number (std::size_t len);
inline cl_SV_number::cl_SV_number (std::size_t len)
	: cl_SV<cl_number,cl_SV_any> (cl_make_heap_SV_number(len)) {}

// Private pointer manipulations. Never throw away a `struct cl_heap_SV_number *'!
inline cl_SV_number::operator cl_heap_SV_number* () const
{
	cl_heap_SV_number* hpointer = (cl_heap_SV_number*)pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
extern const cl_SV_number cl_null_SV_number;
inline cl_SV_number::cl_SV_number ()
	: cl_SV<cl_number,cl_SV_any> ((cl_heap_SV_number*) cl_null_SV_number) {}
class cl_SV_number_init_helper
{
	static int count;
public:
	cl_SV_number_init_helper();
	~cl_SV_number_init_helper();
};
static cl_SV_number_init_helper cl_SV_number_init_helper_instance;

// Copy a simple vector.
inline const cl_SV_number copy (const cl_SV_number& vector)
{ return The(cl_SV_number) (copy((const cl_SV_any&) vector)); }

// Debugging support.
#ifdef CL_DEBUG
extern int cl_SV_number_debug_module;
CL_FORCE_LINK(cl_SV_number_debug_dummy, cl_SV_number_debug_module)
#endif

}  // namespace cln

#endif /* _CL_SV_NUMBER_H */
