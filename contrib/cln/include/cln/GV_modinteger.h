// General vectors of modular integers.

#ifndef _CL_GV_MODINTEGER_H
#define _CL_GV_MODINTEGER_H

#include "cln/modinteger.h"
#include "cln/GV_integer.h"

namespace cln {

// A vector of modular integers (over the same modular integer ring)
// is just a normal vector of integers, with maxbits() operation.

template <>
struct cl_heap_GV<_cl_MI> : cl_heap {
	cl_GV_inner<_cl_MI> v;
	// here room for the elements
};
typedef cl_heap_GV<_cl_MI> cl_heap_GV_MI;

struct cl_GV_MI : public cl_GV<_cl_MI,cl_GV_any> {
public:
	// Constructors.
	cl_GV_MI ();
	cl_GV_MI (const cl_GV_MI&);
	// Create a vector of modular integers.
	cl_GV_MI (std::size_t len, cl_heap_modint_ring* R);
	// Assignment operators.
	cl_GV_MI& operator= (const cl_GV_MI&);
	// Number m of bits allowed per element (-1 if unconstrained).
	sintC maxbits () const
	{
		return ((const cl_heap_GV_I *) pointer)->maxbits();
	}
};
inline cl_GV_MI::cl_GV_MI (const cl_GV_MI& x) : cl_GV<_cl_MI,cl_GV_any> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_GV_MI,cl_GV_MI)
inline cl_GV_MI::cl_GV_MI ()
	: cl_GV<_cl_MI,cl_GV_any> ((cl_heap_GV_MI*) (cl_heap_GV_I*) cl_null_GV_I) {}
inline cl_GV_MI::cl_GV_MI (std::size_t len, cl_heap_modint_ring* R)
	: cl_GV<_cl_MI,cl_GV_any> ((cl_heap_GV_MI*) cl_make_heap_GV_I(len,R->bits)) {}

// Copy a vector.
inline const cl_GV_MI copy (const cl_GV_MI& vector)
{
	return The(cl_GV_MI) (copy((const cl_GV_I&) vector));
}

}  // namespace cln

#endif /* _CL_GV_MODINTEGER_H */
