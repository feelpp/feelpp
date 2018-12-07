// cl_R equal_hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "base/cl_N.h"
#include "real/cl_R.h"
#include "rational/cl_RA.h"
#include "integer/cl_I.h"
#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/misc/cl_SF_eqhashcode.cc"
#include "float/ffloat/misc/cl_FF_eqhashcode.cc"
#include "float/dfloat/misc/cl_DF_eqhashcode.cc"
#include "float/lfloat/misc/cl_LF_eqhashcode.cc"
#include "base/cl_inline2.h"
#include "rational/misc/cl_RA_eqhashcode.cc"

namespace cln {

uint32 CL_FLATTEN equal_hashcode (const cl_R& x)
GEN_R_OP1_7(x, equal_hashcode_inline, return)

}  // namespace cln
