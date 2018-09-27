// minusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#define minusp minusp_inline
#include "rational/cl_RA.h"
#include "integer/cl_I.h"
#undef minusp
#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_minusp.cc"
#include "float/ffloat/elem/cl_FF_minusp.cc"
#include "float/dfloat/elem/cl_DF_minusp.cc"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

bool CL_FLATTEN minusp (const cl_R& x)
#if 0
GEN_R_OP1_2(x, minusp, return)
#else // fully inlined, faster
GEN_R_OP1_7(x, minusp_inline, return)
#endif

}  // namespace cln
