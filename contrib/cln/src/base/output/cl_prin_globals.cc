// Global variables in CLN

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/output.h"


// Implementation.

namespace cln {

cl_print_flags default_print_flags;

int cl_prin_globals_init_helper::count = 0;

cl_prin_globals_init_helper::cl_prin_globals_init_helper()
{
	if (count++ == 0)
		new ((void *)&default_print_flags) cl_print_flags();
}

cl_prin_globals_init_helper::~cl_prin_globals_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up.
	}
}


#if 0 // The default constructors already do this.
AT_INITIALIZATION(default_print_flags)
{
	default_print_flags.rational_base = 10;
	default_print_flags.rational_readably = false;
	default_print_flags.float_readably = false;
	default_print_flags.default_float_format = float_format_ffloat;
	default_print_flags.complex_readably = false;
	default_print_flags.vector_syntax = vsyntax_pretty;
	default_print_flags.univpoly_varname = "x";
}
#endif

}  // namespace cln

