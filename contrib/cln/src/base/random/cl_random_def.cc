// default_random_state.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/random.h"


// Implementation.

namespace cln {
	
random_state default_random_state;

int cl_random_def_init_helper::count = 0;
cl_random_def_init_helper::cl_random_def_init_helper()
{
	if (count++ == 0) {
		default_random_state = random_state();
	}
}

cl_random_def_init_helper::~cl_random_def_init_helper()
{
	if (--count == 0) {
		// Nothing to clean up?
	}
}

}  // namespace cln

