// cl_SV_number debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/output.h"
#include "cln/SV_number.h"
#include "cln/io.h"
#include "vector/cl_SV_io.h"

namespace cln {

static void print_for_debug (std::ostream& stream, const cl_print_flags& flags, const cl_number& z)
{
	unused stream; // must be cl_debugout
	unused flags; // must be default_print_flags
	z.debug_print();
}

static void dprint (cl_heap* pointer)
{
	var const cl_SV_number& obj = *(const cl_SV_number*)&pointer;
	fprint(cl_debugout, "(cl_SV_number) ");
	print_vector(cl_debugout,default_print_flags,&print_for_debug,obj);
}
AT_INITIALIZATION(dprint_SV_number)
{ extern cl_class& cl_class_svector_number(); cl_class_svector_number().dprint = dprint; }

// This dummy links in this module when <cln/SV_number.h> requires it.
int cl_SV_number_debug_module;

}  // namespace cln
