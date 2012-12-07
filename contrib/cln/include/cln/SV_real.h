// Simple vectors of real numbers.

#ifndef _CL_SV_REAL_H
#define _CL_SV_REAL_H

#include "cln/number.h"
#include "cln/SV_complex.h"
#include "cln/real_class.h"
#include "cln/io.h"

namespace cln {

// A vector of real numbers is just a normal vector of numbers.

typedef cl_heap_SV<cl_R> cl_heap_SV_R;

struct cl_SV_R : public cl_SV<cl_R,cl_SV_N> {
public:
	// Constructors.
	cl_SV_R () : cl_SV<cl_R,cl_SV_N> ((cl_heap_SV_R*) (cl_heap_SV_number*) cl_null_SV_number) {};
	cl_SV_R (const cl_SV_R&);
	explicit cl_SV_R (std::size_t len) : cl_SV<cl_R,cl_SV_N> ((cl_heap_SV_R*) cl_make_heap_SV_number(len)) {};
	// Assignment operators.
	cl_SV_R& operator= (const cl_SV_R&);
	// Private pointer manipulations.
	cl_SV_R (cl_heap_SV_R* p) : cl_SV<cl_R,cl_SV_N> (p) {}
	cl_SV_R (cl_private_thing p) : cl_SV<cl_R,cl_SV_N> (p) {}
};
inline cl_SV_R::cl_SV_R (const cl_SV_R& x) : cl_SV<cl_R,cl_SV_N> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SV_R,cl_SV_R)

// Copy a simple vector.
inline const cl_SV_R copy (const cl_SV_R& vector)
{
	return The(cl_SV_R) (copy((const cl_SV_N&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_SV_R& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_SV_number& vector);
	extern void print_real (std::ostream& stream, const cl_print_flags& flags, const cl_R& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_R&))
	             &print_real,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_SV_R)

}  // namespace cln

#endif /* _CL_SV_REAL_H */
