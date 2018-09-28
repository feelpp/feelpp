// cl_class_ratio.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational.h"


// Implementation.

#include "rational/cl_RA.h"

namespace cln {

static void ratio_destructor (cl_heap* pointer)
{
	(*(cl_heap_ratio*)pointer).~cl_heap_ratio();
}

cl_class cl_class_ratio = {
	ratio_destructor,
	cl_class_flags_subclass_complex | cl_class_flags_subclass_real | cl_class_flags_subclass_rational
};

}  // namespace cln
