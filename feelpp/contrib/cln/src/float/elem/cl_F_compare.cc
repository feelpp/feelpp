// compare().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"
#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

cl_signean compare (const cl_F& x, const cl_F& y)
GEN_F_OP2(x,y, compare, 0, 1, return)

}  // namespace cln
