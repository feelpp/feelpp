// round1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "cln/float.h"

namespace cln {

const cl_I round1 (const cl_R& x)
GEN_R_OP1_2(x, round1, return)

}  // namespace cln
