// cl_FF debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/ffloat.h"
#include "cln/io.h"
#include "cln/float_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_FF& obj = *(const cl_FF*)&pointer;
	fprint(cl_debugout, "(cl_FF) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_FF)
{ cl_register_type_printer(cl_class_ffloat,dprint); }

// This dummy links in this module when <cln/ffloat.h> requires it.
int cl_FF_debug_module;

}  // namespace cln
