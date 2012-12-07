// cl_SF_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

inline bool cl_SF_p (const cl_number& x)
{
	if (!x.pointer_p())
		if (cl_tag((x).word) == cl_SF_tag)
			return true;
	return false;
}

const cl_SF& cl_SF_As (const cl_number& x, const char * filename, int line)
{
	if (cl_SF_p(x)) {
		DeclareType(cl_SF,x);
		return x;
	} else
		throw as_exception(x,"a short-float number",filename,line);
}

}  // namespace cln
