// cl_C debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/complex.h"
#include "cln/io.h"
#include "cln/complex_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_N& obj = *(const cl_N*)&pointer;
	fprint(cl_debugout, "(cl_N) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_N)
{ cl_register_type_printer(cl_class_complex,dprint); }

// This dummy links in this module when <cln/complex.h> requires it.
int cl_C_debug_module;

extern int cl_R_debug_module;
static void* dummy[] = { &dummy,
	&cl_R_debug_module
};

}  // namespace cln
