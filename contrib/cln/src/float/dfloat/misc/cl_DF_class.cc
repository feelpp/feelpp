// cl_class_dfloat.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

namespace cln {

cl_class cl_class_dfloat = {
	NULL,		// empty destructor
	cl_class_flags_subclass_complex | cl_class_flags_subclass_real | cl_class_flags_subclass_float
};

}  // namespace cln
