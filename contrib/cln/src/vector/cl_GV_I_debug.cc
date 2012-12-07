// cl_GV_I debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/GV_integer.h"
#include "cln/io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_GV_I& obj = *(const cl_GV_I*)&pointer;
	fprint(cl_debugout, "(cl_GV_I) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_GV_I)
{ extern cl_class& cl_class_gvector_integer(); cl_class_gvector_integer().dprint = dprint; }

// This dummy links in this module when <cln/GV_integer.h> requires it.
int cl_GV_I_debug_module;

}  // namespace cln
