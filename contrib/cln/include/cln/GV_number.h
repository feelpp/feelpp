// General vectors of numbers.

#ifndef _CL_GV_NUMBER_H
#define _CL_GV_NUMBER_H

#include "cln/number.h"
#include "cln/GV.h"

namespace cln {

typedef cl_heap_GV<cl_number> cl_heap_GV_number;

struct cl_GV_number : public cl_GV<cl_number,cl_GV_any> {
public:
	// Constructors.
	cl_GV_number ();
	cl_GV_number (const cl_GV_number&);
	explicit cl_GV_number (std::size_t len);
	// Assignment operators.
	cl_GV_number& operator= (const cl_GV_number&);
	// Private pointer manipulations.
	cl_GV_number (cl_heap_GV_number* p) : cl_GV<cl_number,cl_GV_any> (p) {}
	cl_GV_number (cl_private_thing p) : cl_GV<cl_number,cl_GV_any> (p) {}
};
inline cl_GV_number::cl_GV_number (const cl_GV_number& x) : cl_GV<cl_number,cl_GV_any> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_GV_number,cl_GV_number)
extern cl_heap_GV_number* cl_make_heap_GV_number (std::size_t len);
inline cl_GV_number::cl_GV_number (std::size_t len)
	: cl_GV<cl_number,cl_GV_any> (cl_make_heap_GV_number(len)) {}

// Private pointer manipulations. Never throw away a `struct cl_heap_GV_number *'!
extern const cl_GV_number cl_null_GV_number;
inline cl_GV_number::cl_GV_number ()
	: cl_GV<cl_number,cl_GV_any> ((cl_heap_GV_number*) cl_null_GV_number) {}
class cl_GV_number_init_helper
{
	static int count;
public:
	cl_GV_number_init_helper();
	~cl_GV_number_init_helper();
};
static cl_GV_number_init_helper cl_GV_number_init_helper_instance;

// Copy a vector.
extern const cl_GV_number copy (const cl_GV_number&);

// Debugging support.
#ifdef CL_DEBUG
extern int cl_GV_number_debug_module;
CL_FORCE_LINK(cl_GV_number_debug_dummy, cl_GV_number_debug_module)
#endif

}  // namespace cln

#endif /* _CL_GV_NUMBER_H */
