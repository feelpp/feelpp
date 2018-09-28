// recip().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"

namespace cln {

ALL_cl_LF_OPERATIONS_SAME_PRECISION()

const cl_LF recip (const cl_LF& x)
{
	return encode_LF1(TheLfloat(x)->len) / x;
}

}  // namespace cln
