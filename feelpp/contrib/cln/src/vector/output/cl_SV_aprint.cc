// fprint() for cl_SV_ringelt.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/SV_ringelt.h"


// Implementation.

#include "cln/output.h"

namespace cln {

void fprint (std::ostream& stream, const cl_ring& R, const cl_SV_ringelt& vector)
{
	const cl_print_flags& flags = default_print_flags;
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
		R->_fprint(stream,vector[i]);
	}
	if (flags.vector_syntax == vsyntax_commonlisp)
		fprintchar(stream,')');
	else
		fprintchar(stream,']');
}

}  // namespace cln
