// cl_I_As().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "base/cl_N.h"

namespace cln {

// Cf. cl_I_p in cl_I_ring.cc.
// But here, for better inlining in g++, it is preferrable to finish every
// alternative with either "return true;" or "return false;".

inline bool cl_I_p (const cl_number& x)
{
	if (!x.pointer_p())
		switch (x.nonpointer_tag()) {
		case cl_FN_tag:
			return true;
		}
	else
		if (x.pointer_type() == &cl_class_bignum)
			return true;
	return false;
}

const cl_I& cl_I_As (const cl_number& x, const char * filename, int line)
{
	if (cl_I_p(x)) {
		DeclareType(cl_I,x);
		return x;
	} else
		throw as_exception(x,"an integer",filename,line);
}

}  // namespace cln
