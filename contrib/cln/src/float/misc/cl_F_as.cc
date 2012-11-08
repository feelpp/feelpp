// cl_F_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

inline bool cl_F_p (const cl_number& x)
{
	if (!x.pointer_p())
		switch (cl_tag((x).word)) {
		case cl_SF_tag:
		#if defined(CL_WIDE_POINTERS)
		case cl_FF_tag:
		#endif
			return true;
		}
	else
		if (x.heappointer->type->flags & cl_class_flags_subclass_float)
			return true;
	return false;
}

const cl_F& cl_F_As (const cl_number& x, const char * filename, int line)
{
	if (cl_F_p(x)) {
		DeclareType(cl_F,x);
		return x;
	} else
		throw as_exception(x,"a floating-point number",filename,line);
}

}  // namespace cln
