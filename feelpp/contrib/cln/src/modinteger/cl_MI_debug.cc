// cl_modint_ring debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/modinteger.h"
#include "cln/io.h"
#include "cln/integer_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_modint_ring& obj = *(const cl_modint_ring*)&pointer;
	fprint(cl_debugout, "(cl_modint_ring) Z mod ");
	fprint(cl_debugout, obj->modulus);
}
AT_INITIALIZATION(dprint_modint_ring)
{ cl_register_type_printer(cl_class_modint_ring,dprint); }

void cl_MI::debug_print () const
{
	fprint(cl_debugout, *this);
	fprint(cl_debugout, "\n");
}

// This dummy links in this module when <cln/modinteger.h> requires it.
int cl_MI_debug_module;

}  // namespace cln
