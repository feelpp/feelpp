// cl_class_complex.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/complex.h"


// Implementation.

#include "complex/cl_C.h"

namespace cln {

static void complex_destructor (cl_heap* pointer)
{
	(*(cl_heap_complex*)pointer).~cl_heap_complex();
}

cl_class cl_class_complex = {
	complex_destructor,
	cl_class_flags_subclass_complex
};

}  // namespace cln
