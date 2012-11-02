// cl_SF debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/sfloat.h"
#include "cln/io.h"
#include "cln/float_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_SF& obj = *(const cl_SF*)&pointer;
	fprint(cl_debugout, "(cl_SF) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_SF)
{ cl_register_type_printer(cl_class_sfloat,dprint); }

// This dummy links in this module when <cln/sfloat.h> requires it.
int cl_SF_debug_module;

}  // namespace cln
