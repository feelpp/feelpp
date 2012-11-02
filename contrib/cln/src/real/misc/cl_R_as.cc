// cl_R_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

// Cf. cl_R_p in cl_R_ring.cc.
// But here, for better inlining in g++, it is preferrable to finish every
// alternative with either "return true;" or "return false;".

inline bool cl_R_p (const cl_number& x)
{
	if (!x.pointer_p())
		switch (x.nonpointer_tag()) {
		case cl_FN_tag:
		case cl_SF_tag:
		#if defined(CL_WIDE_POINTERS)
		case cl_FF_tag:
		#endif
			return true;
		}
	else
		if (x.pointer_type()->flags & cl_class_flags_subclass_real)
			return true;
	return false;
}

const cl_R& cl_R_As (const cl_number& x, const char * filename, int line)
{
	if (cl_R_p(x)) {
		DeclareType(cl_R,x);
		return x;
	} else
		throw as_exception(x,"a real number",filename,line);
}

}  // namespace cln
