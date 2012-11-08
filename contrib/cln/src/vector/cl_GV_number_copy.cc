// copy().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#define CL_GV_NO_RANGECHECKS
#include "cln/GV_number.h"


// Implementation.

namespace cln {

const cl_GV_number copy (const cl_GV_number& v)
{
	std::size_t len = v.size();
	cl_GV_number w = cl_GV_number(len);
	cl_GV_number::copy_elements(v, 0, w, 0, len);
	return w;
}

}  // namespace cln
