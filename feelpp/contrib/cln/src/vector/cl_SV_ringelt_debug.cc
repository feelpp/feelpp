// cl_SV_ringelt debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/SV_ringelt.h"
#include "cln/io.h"

namespace cln {

extern void cl_dprint_unknown (cl_heap* pointer);
static void dprint (cl_heap* pointer)
{
//	var const cl_SV_ringelt& obj = *(const cl_SV_ringelt*)&pointer;
	fprint(cl_debugout, "(cl_SV_ringelt) ");
	cl_dprint_unknown(pointer);
}
AT_INITIALIZATION(dprint_SV_ringelt)
{ extern cl_class& cl_class_svector_ringelt(); cl_class_svector_ringelt().dprint = dprint; }

// This dummy links in this module when <cln/SV_ringelt.h> requires it.
int cl_SV_ringelt_debug_module;

}  // namespace cln
