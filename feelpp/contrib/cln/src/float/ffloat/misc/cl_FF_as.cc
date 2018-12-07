// cl_FF_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/ffloat.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

inline bool cl_FF_p (const cl_number& x)
{
#if defined(CL_WIDE_POINTERS)
	if (!x.pointer_p())
		if (cl_tag((x).word) == cl_FF_tag)
			return true;
#else
	if (x.pointer_p())
		if (x.heappointer->type == &cl_class_ffloat)
			return true;
#endif
	return false;
}

const cl_FF& cl_FF_As (const cl_number& x, const char * filename, int line)
{
	if (cl_FF_p(x)) {
		DeclareType(cl_FF,x);
		return x;
	} else
		throw as_exception(x,"a single-float number",filename,line);
}

}  // namespace cln
