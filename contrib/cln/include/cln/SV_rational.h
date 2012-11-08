// Simple vectors of rational numbers.

#ifndef _CL_SV_RATIONAL_H
#define _CL_SV_RATIONAL_H

#include "cln/number.h"
#include "cln/SV_real.h"
#include "cln/rational_class.h"
#include "cln/io.h"

namespace cln {

// A vector of rational numbers is just a normal vector of real numbers.

typedef cl_heap_SV<cl_RA> cl_heap_SV_RA;

struct cl_SV_RA : public cl_SV<cl_RA,cl_SV_R> {
public:
	// Constructors.
	cl_SV_RA () : cl_SV<cl_RA,cl_SV_R> ((cl_heap_SV_RA*) (cl_heap_SV_number*) cl_null_SV_number) {};
	cl_SV_RA (const cl_SV_RA&);
	explicit cl_SV_RA (std::size_t len) : cl_SV<cl_RA,cl_SV_R> ((cl_heap_SV_RA*) cl_make_heap_SV_number(len)) {};
	// Assignment operators.
	cl_SV_RA& operator= (const cl_SV_RA&);
	// Private pointer manipulations.
	cl_SV_RA (cl_heap_SV_RA* p) : cl_SV<cl_RA,cl_SV_R> (p) {}
	cl_SV_RA (cl_private_thing p) : cl_SV<cl_RA,cl_SV_R> (p) {}
};
inline cl_SV_RA::cl_SV_RA (const cl_SV_RA& x) : cl_SV<cl_RA,cl_SV_R> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SV_RA,cl_SV_RA)

// Copy a simple vector.
inline const cl_SV_RA copy (const cl_SV_RA& vector)
{
	return The(cl_SV_RA) (copy((const cl_SV_R&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_SV_RA& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_SV_number& vector);
	extern void print_rational (std::ostream& stream, const cl_print_flags& flags, const cl_RA& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_RA&))
	             &print_rational,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_SV_RA)

}  // namespace cln

#endif /* _CL_SV_RAATIONAL_H */
