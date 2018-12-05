// General vectors of rational numbers.

#ifndef _CL_GV_RATIONAL_H
#define _CL_GV_RATIONAL_H

#include "cln/number.h"
#include "cln/GV_real.h"
#include "cln/rational_class.h"
#include "cln/io.h"

namespace cln {

// A vector of rational numbers is just a normal vector of real numbers.

typedef cl_heap_GV<cl_RA> cl_heap_GV_RA;

struct cl_GV_RA : public cl_GV<cl_RA,cl_GV_R> {
public:
	// Constructors.
	cl_GV_RA ();
	cl_GV_RA (const cl_GV_RA&);
	explicit cl_GV_RA (std::size_t len);
	// Assignment operators.
	cl_GV_RA& operator= (const cl_GV_RA&);
	// Private pointer manipulations.
	cl_GV_RA (cl_heap_GV_RA* p) : cl_GV<cl_RA,cl_GV_R> (p) {}
	cl_GV_RA (cl_private_thing p) : cl_GV<cl_RA,cl_GV_R> (p) {}
};
inline cl_GV_RA::cl_GV_RA (const cl_GV_RA& x) : cl_GV<cl_RA,cl_GV_R> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_GV_RA,cl_GV_RA)
inline cl_GV_RA::cl_GV_RA (std::size_t len)
	: cl_GV<cl_RA,cl_GV_R> ((cl_heap_GV_RA*) cl_make_heap_GV_number(len)) {}
inline cl_GV_RA::cl_GV_RA ()
	: cl_GV<cl_RA,cl_GV_R> ((cl_heap_GV_RA*) (cl_heap_GV_number*) cl_null_GV_number) {}

// Copy a vector.
inline const cl_GV_RA copy (const cl_GV_RA& vector)
{
	return The(cl_GV_RA) (copy((const cl_GV_R&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_GV_RA& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_GV_number& vector);
	extern void print_rational (std::ostream& stream, const cl_print_flags& flags, const cl_RA& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_RA&))
	             &print_rational,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_GV_RA)

}  // namespace cln

#endif /* _CL_GV_RAATIONAL_H */
