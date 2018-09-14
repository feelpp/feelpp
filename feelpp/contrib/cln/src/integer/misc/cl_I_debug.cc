// cl_I debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

#include "cln/integer.h"
#include "cln/io.h"
#include "cln/integer_io.h"

namespace cln {

static void dprint (cl_heap* pointer)
{
	var const cl_I& obj = *(const cl_I*)&pointer;
	fprint(cl_debugout, "(cl_I) ");
	fprint(cl_debugout, obj);
}
AT_INITIALIZATION(dprint_I)
{
	cl_register_type_printer(cl_class_fixnum,dprint);
	cl_register_type_printer(cl_class_bignum,dprint);
}

// This dummy links in this module when <cln/integer.h> requires it.
int cl_I_debug_module;

}  // namespace cln
