// cl_univpoly_ring debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/univpoly.h"
#include "polynomial/cl_UP.h"
#include "cln/io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
//	var const cl_univpoly_ring& obj = *(const cl_univpoly_ring*)&pointer;
	fprint(cl_debugout, "(cl_univpoly_ring) ring");
	fprint(cl_debugout, get_varname((cl_heap_univpoly_ring*)pointer));
}
AT_INITIALIZATION(dprint_univpoly_ring)
{ cl_register_type_printer(cl_class_univpoly_ring,dprint); }

void cl_UP::debug_print () const
{
	fprint(cl_debugout, *this);
	fprint(cl_debugout, "\n");
}

// This dummy links in this module when <cl_univpoly_ring.h> requires it.
int cl_UP_debug_module;

}  // namespace cln
