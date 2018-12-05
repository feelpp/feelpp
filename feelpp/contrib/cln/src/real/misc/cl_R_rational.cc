// rational().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"

namespace cln {

const cl_RA rational (const cl_R& x)
GEN_R_OP1_2(x, rational, return)

}  // namespace cln
