// cl_class_fixnum.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

namespace cln {

cl_class cl_class_fixnum = {
	NULL,		// destructor not used, since not heap objects
	cl_class_flags_subclass_complex | cl_class_flags_subclass_real | cl_class_flags_subclass_rational
};

AT_INITIALIZATION(ini_class_fixnum)
{
	cl_immediate_classes[cl_FN_tag] = &cl_class_fixnum;
}

}  // namespace cln
