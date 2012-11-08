// Input functions.

#ifndef _CL_INPUT_H
#define _CL_INPUT_H

#include "cln/types.h"
#include "cln/floatformat.h"
#include "cln/io.h"

namespace cln {

struct cl_read_float_flags {
	// The float format used when reading floats with exponent marker 'E'.
	float_format_t default_float_format;
	// The float format used when reading floats with exponent marker 'L'.
	float_format_t default_lfloat_format;
	// Flag whether floats specified with more digits than corresponding
	// to the exponent marker they contain, but without _nnn suffix, will
	// get a precision corresponding to their number of significant digits.
	bool mantissa_dependent_float_format;
};

// Specifies the possible results of a read operation.
enum cl_read_syntax_t {
	syntax_integer = 1 << 0,				// -> cl_I
	syntax_ratio = 1 << 1,					// -> cl_RA
	syntax_rational = syntax_integer | syntax_ratio,	// -> cl_RA
	syntax_sfloat = 1 << 2,					// -> cl_SF
	syntax_ffloat = 1 << 3,					// -> cl_FF
	syntax_dfloat = 1 << 4,					// -> cl_DF
	syntax_lfloat = 1 << 5,					// -> cl_LF
	syntax_float = syntax_sfloat | syntax_ffloat | syntax_dfloat | syntax_lfloat,
								// -> cl_F
	syntax_real = syntax_rational | syntax_float,		// -> cl_R
	syntax_complex = 1 << 6,				// -> cl_N
	syntax_number = syntax_real | syntax_complex,		// -> cl_N
	syntax_maybe_bad = 1 << 7				// avoid errors
};

// Specifies the syntax to be applied to a read operation.
enum cl_read_lsyntax_t {
		// Standard algebraic notation.
	lsyntax_standard = 0,
		// Extended algebraic notation: x+yi
	lsyntax_algebraic = 1 << 0,
		// Common Lisp notation: #b, #o, #x, #r, #c
	lsyntax_commonlisp = 1 << 1,
		// All of them.
	lsyntax_all = lsyntax_algebraic | lsyntax_commonlisp
};

struct cl_read_flags {
	cl_read_syntax_t syntax;
	cl_read_lsyntax_t lsyntax;
	unsigned int rational_base;
	cl_read_float_flags float_flags;
};

}  // namespace cln

#endif /* _CL_INPUT_H */
