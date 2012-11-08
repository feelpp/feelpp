// signum().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"

/* use inline versions of zerop */
#include "base/cl_inline.h"
#include "integer/misc/cl_I_signum.cc"
#include "rational/misc/cl_RA_signum.cc"
/* use inline versions of signum */
#include "base/cl_inline2.h"
#include "float/sfloat/misc/cl_SF_signum.cc"
#include "float/ffloat/misc/cl_FF_signum.cc"
#include "float/dfloat/misc/cl_DF_signum.cc"
#include "float/lfloat/misc/cl_LF_signum.cc"

namespace cln {

const cl_R CL_FLATTEN signum (const cl_R& x)
GEN_R_OP1_7(x, signum_inline, return)

}  // namespace cln
