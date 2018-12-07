// print_vector().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/GV_complex.h"
#include "cln/GV_real.h"
#include "cln/GV_rational.h"
#include "cln/GV_integer.h"
#include "vector/cl_GV_io.h"


// Implementation.

#include "cln/output.h"

namespace cln {

void print_vector (std::ostream& stream, const cl_print_flags& flags, void (* printfun) (std::ostream&, const cl_print_flags&, const cl_number&), const cl_GV_number& vector)
{
	std::size_t len = vector.size();
	if (flags.vector_syntax == vsyntax_commonlisp) {
		fprintchar(stream,'#');
		fprintchar(stream,'(');
	} else
		fprintchar(stream,'[');
	for (std::size_t i = 0; i < len; i++) {
		if (i > 0) {
			if (flags.vector_syntax == vsyntax_algebraic)
				fprintchar(stream,',');
			fprintchar(stream,' ');
		}
		// The conversion to cl_number below is needed for SGI CC.
		printfun(stream,flags,(cl_number)vector[i]);
	}
	if (flags.vector_syntax == vsyntax_commonlisp)
		fprintchar(stream,')');
	else
		fprintchar(stream,']');
}

}  // namespace cln
