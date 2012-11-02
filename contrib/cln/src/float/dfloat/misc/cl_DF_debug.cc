// cl_DF debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/dfloat.h"
#include "cln/io.h"
#include "cln/float_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_DF& obj = *(const cl_DF*)&pointer;
	fprint(cl_debugout, "(cl_DF) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_DF)
{ cl_register_type_printer(cl_class_dfloat,dprint); }

// This dummy links in this module when <cln/dfloat.h> requires it.
int cl_DF_debug_module;

}  // namespace cln
