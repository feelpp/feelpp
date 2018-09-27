// General vectors of complex numbers.

#ifndef _CL_GV_COMPLEX_H
#define _CL_GV_COMPLEX_H

#include "cln/number.h"
#include "cln/GV_number.h"
#include "cln/complex_class.h"
#include "cln/io.h"

namespace cln {

// A vector of complex numbers is just a normal vector of numbers.

typedef cl_heap_GV<cl_N> cl_heap_GV_N;

struct cl_GV_N : public cl_GV<cl_N,cl_GV_number> {
public:
	// Constructors.
	cl_GV_N ();
	cl_GV_N (const cl_GV_N&);
	explicit cl_GV_N (std::size_t len);
	// Assignment operators.
	cl_GV_N& operator= (const cl_GV_N&);
	// Private pointer manipulations.
	cl_GV_N (cl_heap_GV_N* p) : cl_GV<cl_N,cl_GV_number> (p) {}
	cl_GV_N (cl_private_thing p) : cl_GV<cl_N,cl_GV_number> (p) {}
};
inline cl_GV_N::cl_GV_N (const cl_GV_N& x) : cl_GV<cl_N,cl_GV_number> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_GV_N,cl_GV_N)
inline cl_GV_N::cl_GV_N (std::size_t len)
	: cl_GV<cl_N,cl_GV_number> ((cl_heap_GV_N*) cl_make_heap_GV_number(len)) {}
inline cl_GV_N::cl_GV_N ()
	: cl_GV<cl_N,cl_GV_number> ((cl_heap_GV_N*) (cl_heap_GV_number*) cl_null_GV_number) {}

// Copy a vector.
inline const cl_GV_N copy (const cl_GV_N& vector)
{
	return The(cl_GV_N) (copy((const cl_GV_number&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_GV_N& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_GV_number& vector);
	extern void print_complex (std::ostream& stream, const cl_print_flags& flags, const cl_N& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_N&))
	             &print_complex,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_GV_N)

}  // namespace cln

#endif /* _CL_GV_COMPLEX_H */
