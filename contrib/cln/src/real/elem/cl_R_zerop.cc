// zerop().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#define zerop zerop_inline
#include "real/cl_R.h"
#include "rational/cl_RA.h"
#include "integer/cl_I.h"
#undef zerop
#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_zerop.cc"
#include "float/ffloat/elem/cl_FF_zerop.cc"
#include "float/dfloat/elem/cl_DF_zerop.cc"
#include "float/lfloat/elem/cl_LF_zerop.cc"

namespace cln {

bool CL_FLATTEN zerop (const cl_R& x)
#if 0
GEN_R_OP1_2(x, zerop, return)
#else // fully inlined, faster
GEN_R_OP1_7(x, zerop_inline, return)
#endif

}  // namespace cln
