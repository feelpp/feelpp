// cl_ring debugging support.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ring.h"


// Implementation.

#include "cln/io.h"

namespace cln {

void cl_ring_element::debug_print () const
{
	fprint(cl_debugout, *this);
	cl_debugout << std::endl; // newline and flush output
}

// This dummy links in this module when <cln/ring.h> requires it.
int cl_ring_debug_module;

}  // namespace cln
