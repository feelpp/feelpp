// division_by_0_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/exception.h"

namespace cln {

division_by_0_exception::division_by_0_exception ()
	: runtime_exception("Division by zero.")
{}

}  // namespace cln
