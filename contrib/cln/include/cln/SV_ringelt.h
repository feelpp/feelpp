// Simple vectors of ring elements.

#ifndef _CL_SV_RINGELT_H
#define _CL_SV_RINGELT_H

#include "cln/ring.h"
#include "cln/SV.h"
#include "cln/io.h"

namespace cln {

typedef cl_heap_SV<_cl_ring_element> cl_heap_SV_ringelt;

struct cl_SV_ringelt : public cl_SV<_cl_ring_element,cl_SV_any> {
public:
	// Constructors.
	cl_SV_ringelt ();
	cl_SV_ringelt (const cl_SV_ringelt&);
	explicit cl_SV_ringelt (std::size_t len);
	// Assignment operators.
	cl_SV_ringelt& operator= (const cl_SV_ringelt&);
	// Private pointer manipulations.
	operator cl_heap_SV_ringelt* () const;
	cl_SV_ringelt (cl_heap_SV_ringelt* p) : cl_SV<_cl_ring_element,cl_SV_any> (p) {}
	cl_SV_ringelt (cl_private_thing p) : cl_SV<_cl_ring_element,cl_SV_any> (p) {}
};
inline cl_SV_ringelt::cl_SV_ringelt (const cl_SV_ringelt& x) : cl_SV<_cl_ring_element,cl_SV_any> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SV_ringelt,cl_SV_ringelt)
// Returns a new simple vector with uninitialized contents.
extern cl_heap_SV_ringelt* cl_make_heap_SV_ringelt_uninit (std::size_t len);
// Returns a new simple vector with all elements initialized to some value.
extern cl_heap_SV_ringelt* cl_make_heap_SV_ringelt (std::size_t len);
inline cl_SV_ringelt::cl_SV_ringelt (std::size_t len)
	: cl_SV<_cl_ring_element,cl_SV_any> (cl_make_heap_SV_ringelt(len)) {}

// Private pointer manipulations.
// Never throw away a `struct cl_heap_SV_ringelt *'!
inline cl_SV_ringelt::operator cl_heap_SV_ringelt* () const
{
	cl_heap_SV_ringelt* hpointer = (cl_heap_SV_ringelt*)pointer;
	cl_inc_refcount(*this);
	return hpointer;
}
extern const cl_SV_ringelt cl_null_SV_ringelt;
inline cl_SV_ringelt::cl_SV_ringelt ()
	: cl_SV<_cl_ring_element,cl_SV_any> ((cl_heap_SV_ringelt*) cl_null_SV_ringelt) {}

class cl_SV_ringelt_init_helper
{
	static int count;
public:
	cl_SV_ringelt_init_helper();
	~cl_SV_ringelt_init_helper();
};
static cl_SV_ringelt_init_helper cl_SV_ringelt_init_helper_instance;

// Copy a simple vector.
inline const cl_SV_ringelt copy (const cl_SV_ringelt& vector)
{ return The(cl_SV_ringelt) (copy((const cl_SV_any&) vector)); }

// Output.
extern void fprint (std::ostream& stream, const cl_ring& R, const cl_SV_ringelt& x);

// Debugging support.
#ifdef CL_DEBUG
extern int cl_SV_ringelt_debug_module;
CL_FORCE_LINK(cl_SV_ringelt_debug_dummy, cl_SV_ringelt_debug_module)
#endif

}  // namespace cln

#endif /* _CL_SV_RINGELT_H */
