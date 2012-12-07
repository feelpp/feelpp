// destructor ~cl_timing().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/timing.h"


// Implementation.

namespace cln {

cl_timing::~cl_timing ()
{
	report_fn(*this);
}

}  // namespace cln
