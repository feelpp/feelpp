// Simple vectors of integers.

#ifndef _CL_SV_INTEGER_H
#define _CL_SV_INTEGER_H

#include "cln/number.h"
#include "cln/SV_rational.h"
#include "cln/integer_class.h"
#include "cln/io.h"

namespace cln {

// A vector of integers is just a normal vector of rational numbers.

typedef cl_heap_SV<cl_I> cl_heap_SV_I;

struct cl_SV_I : public cl_SV<cl_I,cl_SV_RA> {
public:
	// Constructors.
	cl_SV_I () : cl_SV<cl_I,cl_SV_RA> ((cl_heap_SV_I*) (cl_heap_SV_number*) cl_null_SV_number) {};
	cl_SV_I (const cl_SV_I&);
	explicit cl_SV_I (std::size_t len) : cl_SV<cl_I,cl_SV_RA> ((cl_heap_SV_I*) cl_make_heap_SV_number(len)) {};
	// Assignment operators.
	cl_SV_I& operator= (const cl_SV_I&);
};
inline cl_SV_I::cl_SV_I (const cl_SV_I& x) : cl_SV<cl_I,cl_SV_RA> (as_cl_private_thing(x)) {}
CL_DEFINE_ASSIGNMENT_OPERATOR(cl_SV_I,cl_SV_I)

// Copy a simple vector.
inline const cl_SV_I copy (const cl_SV_I& vector)
{
	return The(cl_SV_I) (copy((const cl_SV_RA&) vector));
}

// Output.
inline void fprint (std::ostream& stream, const cl_SV_I& x)
{
	extern cl_print_flags default_print_flags;
	extern void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* fun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_SV_number& vector);
	extern void print_integer (std::ostream& stream, const cl_print_flags& flags, const cl_I& z);
	print_vector(stream, default_print_flags,
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_number&))
	             (void (*) (std::ostream&, const cl_print_flags&, const cl_I&))
	             &print_integer,
	             x);
}
CL_DEFINE_PRINT_OPERATOR(cl_SV_I)

}  // namespace cln

#endif /* _CL_SV_INTEGER_H */
