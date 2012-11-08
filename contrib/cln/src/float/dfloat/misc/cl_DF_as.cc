// cl_DF_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

inline bool cl_DF_p (const cl_number& x)
{
	if (x.pointer_p())
		if (x.heappointer->type == &cl_class_dfloat)
			return true;
	return false;
}

const cl_DF& cl_DF_As (const cl_number& x, const char * filename, int line)
{
	if (cl_DF_p(x)) {
		DeclareType(cl_DF,x);
		return x;
	} else
		throw as_exception(x,"a double-float number",filename,line);
}

}  // namespace cln
