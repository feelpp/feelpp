// cl_RA debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/rational.h"
#include "cln/io.h"
#include "cln/rational_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_RA& obj = *(const cl_RA*)&pointer;
	fprint(cl_debugout, "(cl_RA) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_RA)
{ cl_register_type_printer(cl_class_ratio,dprint); }

// This dummy links in this module when <cln/rational.h> requires it.
int cl_RA_debug_module;

extern int cl_I_debug_module;
static void* dummy[] = { &dummy,
	&cl_I_debug_module
};

}  // namespace cln
