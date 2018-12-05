// cl_LF debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/lfloat.h"
#include "cln/io.h"
#include "cln/float_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_LF& obj = *(const cl_LF*)&pointer;
	fprint(cl_debugout, "(cl_LF) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_LF)
{ cl_register_type_printer(cl_class_lfloat,dprint); }

// This dummy links in this module when <cln/lfloat.h> requires it.
int cl_LF_debug_module;

}  // namespace cln
