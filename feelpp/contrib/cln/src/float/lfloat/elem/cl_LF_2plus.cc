// binary operator +

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_LF operator+ (const cl_LF& x1, const cl_LF& x2)
{ GEN_LF_OP2(x1,x2,LF_LF_plus_LF,return) }

}  // namespace cln
