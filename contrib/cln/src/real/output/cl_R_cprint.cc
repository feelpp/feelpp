// print_real().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real_io.h"


// Implementation.

#include "cln/real.h"
#include "real/cl_R.h"
#include "cln/rational_io.h"
#include "cln/float_io.h"

namespace cln {

void print_real (std::ostream& stream, const cl_print_real_flags& flags, const cl_R& z)
{
	if (rationalp(z)) {
		DeclareType(cl_RA,z);
		// rationale Zahl
		print_rational(stream,flags,z);
	} else {
		DeclareType(cl_F,z);
		// Float
		print_float(stream,flags,z);
	}
}

}  // namespace cln
