// General vectors of real numbers.

#ifndef _CL_GV_REAL_H
#define _CL_GV_REAL_H

#include "cln/number.h"
#include "cln/GV_complex.h"
#include "cln/real_class.h"
#include "cln/io.h"

namespace cln {

// A vector of real numbers is just a normal vector of numbers.

typedef cl_heap_GV<cl_R> cl_heap_GV_R;

struct cl_GV_R : public cl_GV<cl_R,cl_GV_N> {
public:
	// Constructors.
	cl_GV_R ();
	cl_GV_R (const cl_GV_R&);
	explicit cl_GV_R (std::size_t len);
	// Assignment operators.
	cl_GV_R& operator= (const cl_GV_R&);
	// Private pointer manipulations.
	cl_GV_R (cl_heap_GV_R* p) : cl_GV<cl_R,cl_GV_N> (p) {}
	cl_GV_R (cl_private_thing p) : cl_GV<cl_R,cl_GV_N> (p) {}
};
inline cl_GV_R::cl_GV_R (const cl_GV_R& x) : cl_GV<cl_R,cl_GV_N> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_GV_R,cl_GV_R)
inline cl_GV_R::cl_GV_R (std::size_t len)
	: cl_GV<cl_R,cl_GV_N> ((cl_heap_GV_R*) cl_make_heap_GV_number(len)) {}
inline cl_GV_R::cl_GV_R ()
	: cl_GV<cl_R,cl_GV_N> ((cl_heap_GV_R*) (cl_heap_GV_number*) cl_null_GV_number) {}

// Copy a vector.
inline const cl_GV_R copy (const cl_GV_R& vector)
{
	return The(cl_GV_R) (copy((const cl_GV_N&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_GV_R& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_GV_number& vector);
	extern void print_real (std::ostream& stream, const cl_print_flags& flags, const cl_R& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_R&))
	             &print_real,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_GV_R)

}  // namespace cln

#endif /* _CL_GV_REAL_H */
