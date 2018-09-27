// cl_R debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.


// Implementation.

namespace cln {

// This dummy links in this module when <cln/real.h> requires it.
int cl_R_debug_module;

extern int cl_SF_debug_module;
extern int cl_FF_debug_module;
extern int cl_DF_debug_module;
extern int cl_LF_debug_module;
extern int cl_RA_debug_module;
static void* dummy[] = { &dummy,
	&cl_SF_debug_module,
	&cl_FF_debug_module,
	&cl_DF_debug_module,
	&cl_LF_debug_module,
	&cl_RA_debug_module
};

}  // namespace cln
