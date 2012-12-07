// cl_notify_composite().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "modinteger/cl_MI.h"


// Implementation.

#include "cln/io.h"

namespace cln {

cl_composite_condition* cl_notify_composite (const cl_modint_ring& R, const cl_I& nonunit)
{
	return new cl_composite_condition(R->modulus,gcd(R->modulus,nonunit));
}

}  // namespace cln
