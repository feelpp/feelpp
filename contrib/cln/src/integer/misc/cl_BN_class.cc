// cl_class_bignum.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

namespace cln {

cl_class cl_class_bignum = {
	NULL,		// empty destructor
	cl_class_flags_subclass_complex | cl_class_flags_subclass_real | cl_class_flags_subclass_rational
};

}  // namespace cln
