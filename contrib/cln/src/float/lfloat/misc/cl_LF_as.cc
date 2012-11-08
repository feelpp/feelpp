// cl_LF_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

inline bool cl_LF_p (const cl_number& x)
{
	if (x.pointer_p())
		if (x.heappointer->type == &cl_class_lfloat)
			return true;
	return false;
}

const cl_LF& cl_LF_As (const cl_number& x, const char * filename, int line)
{
	if (cl_LF_p(x)) {
		DeclareType(cl_LF,x);
		return x;
	} else
		throw as_exception(x,"a long-float number",filename,line);
}

}  // namespace cln
