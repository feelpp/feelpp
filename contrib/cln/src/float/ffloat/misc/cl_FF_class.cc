// cl_class_ffloat.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

namespace cln {

cl_class cl_class_ffloat = {
#ifdef CL_WIDE_POINTERS
	NULL,		// destructor not used, since not heap objects
#else
	NULL,		// empty destructor
#endif
	cl_class_flags_subclass_complex | cl_class_flags_subclass_real | cl_class_flags_subclass_float
};

#ifdef CL_WIDE_POINTERS
AT_INITIALIZATION(ini_class_ffloat)
{
	cl_immediate_classes[cl_FF_tag] = &cl_class_ffloat;
}
#endif

}  // namespace cln
