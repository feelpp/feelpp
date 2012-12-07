// Output functions.

#ifndef _CL_OUTPUT_H
#define _CL_OUTPUT_H

#include "cln/types.h"
#include "cln/floatformat.h"
#include "cln/io.h"
#include "cln/string.h"

namespace cln {

struct cl_print_rational_flags {
	// Base in which rational numbers are to be printed.
	unsigned int rational_base;
	// Flag whether to print radix specifiers in Common Lisp syntax for
	// rational numbers (#nR or #b or #o or #x prefixes, trailing dot).
	bool rational_readably;
// Constructor.
	cl_print_rational_flags () :
		rational_base (10),
		rational_readably (false) {}
};

struct cl_print_float_flags {
	// Flag whether to prefer type specific exponent markers over 'E'.
	bool float_readably;
	// If !float_readably, the format which earns the 'E' exponent marker.
	float_format_t default_float_format;
// Constructor.
	cl_print_float_flags () :
		float_readably (false),
		default_float_format (float_format_ffloat) {}
};

struct cl_print_real_flags : cl_print_rational_flags, cl_print_float_flags {};

struct cl_print_complex_flags {
	// Flag whether to use the Common Lisp #C(realpart imagpart) syntax,
	bool complex_readably;
// Constructor.
	cl_print_complex_flags () :
		complex_readably (false) {}
};

struct cl_print_number_flags : cl_print_real_flags, cl_print_complex_flags {};

enum cl_print_vector_syntax_t {
	vsyntax_algebraic,	// [a, b, c]
	vsyntax_pretty,		// [a b c]
	vsyntax_commonlisp	// #(a b c)
};

struct cl_print_vector_flags {
	cl_print_vector_syntax_t vector_syntax;
// Constructor.
	cl_print_vector_flags () :
		vector_syntax (vsyntax_pretty) {}
};

struct cl_print_univpoly_flags {
	cl_string univpoly_varname;
// Constructor.
	cl_print_univpoly_flags () :
		univpoly_varname ("x") {}
};

struct cl_print_flags : cl_print_number_flags, cl_print_vector_flags, cl_print_univpoly_flags {};

extern cl_print_flags default_print_flags;

}  // namespace cln

#endif /* _CL_OUTPUT_H */
