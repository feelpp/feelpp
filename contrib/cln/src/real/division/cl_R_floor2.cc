// floor2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "cln/float.h"
#include "real/division/cl_R_div_t.h"

namespace cln {

const cl_R_div_t floor2 (const cl_R& x)
GEN_R_OP1_2(x, floor2, return)

}  // namespace cln
