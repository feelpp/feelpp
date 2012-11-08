// Simple vectors of complex numbers.

#ifndef _CL_SV_COMPLEX_H
#define _CL_SV_COMPLEX_H

#include "cln/number.h"
#include "cln/SV_number.h"
#include "cln/complex_class.h"
#include "cln/io.h"

namespace cln {

// A vector of complex numbers is just a normal vector of numbers.

typedef cl_heap_SV<cl_N> cl_heap_SV_N;

struct cl_SV_N : public cl_SV<cl_N,cl_SV_number> {
public:
	// Constructors.
	cl_SV_N () : cl_SV<cl_N,cl_SV_number> ((cl_heap_SV_N*) (cl_heap_SV_number*) cl_null_SV_number) {};
	cl_SV_N (const cl_SV_N&);
	explicit cl_SV_N (std::size_t len) : cl_SV<cl_N,cl_SV_number> ((cl_heap_SV_N*) cl_make_heap_SV_number(len)) {};
	// Assignment operators.
	cl_SV_N& operator= (const cl_SV_N&);
	// Private pointer manipulations.
	cl_SV_N (cl_heap_SV_N* p) : cl_SV<cl_N,cl_SV_number> (p) {}
	cl_SV_N (cl_private_thing p) : cl_SV<cl_N,cl_SV_number> (p) {}
};
inline cl_SV_N::cl_SV_N (const cl_SV_N& x) : cl_SV<cl_N,cl_SV_number> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SV_N,cl_SV_N)

// Copy a simple vector.
inline const cl_SV_N copy (const cl_SV_N& vector)
{
	return The(cl_SV_N) (copy((const cl_SV_number&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_SV_N& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_SV_number& vector);
	extern void print_complex (std::ostream& stream, const cl_print_flags& flags, const cl_N& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_N&))
	             &print_complex,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_SV_N)

}  // namespace cln

#endif /* _CL_SV_COMPLEX_H */
