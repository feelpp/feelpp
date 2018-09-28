// cl_GV_number debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/output.h"
#include "cln/GV_number.h"
#include "cln/io.h"
#include "vector/cl_GV_io.h"

namespace cln {

static void print_for_debug (std::ostream& stream, const cl_print_flags& flags, const cl_number& z)
{
	unused stream; // must be cl_debugout
	unused flags; // must be default_print_flags
	z.debug_print();
}

static void dprint (cl_heap* pointer)
{
	var const cl_GV_number& obj = *(const cl_GV_number*)&pointer;
	fprint(cl_debugout, "(cl_GV_number) ");
	print_vector(cl_debugout,default_print_flags,&print_for_debug,obj);
}
AT_INITIALIZATION(dprint_GV_number)
{ extern cl_class& cl_class_gvector_number(); cl_class_gvector_number().dprint = dprint; }

// This dummy links in this module when <cln/GV_number.h> requires it.
int cl_GV_number_debug_module;

extern int cl_GV_I_debug_module;
static void* dummy[] = { &dummy,
	&cl_GV_I_debug_module
};

}  // namespace cln
