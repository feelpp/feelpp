// copy().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#define CL_GV_NO_RANGECHECKS
#include "cln/GV_integer.h"


// Implementation.

namespace cln {

const cl_GV_I copy (const cl_GV_I& v)
{
	std::size_t len = v.size();
	cl_GV_I w = cl_GV_I(len, v.maxbits());
	cl_GV_I::copy_elements(v, 0, w, 0, len);
	return w;
}

}  // namespace cln
