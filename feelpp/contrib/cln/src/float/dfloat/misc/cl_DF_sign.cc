// float_sign().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

#include "base/cl_inline.h"
#include "float/dfloat/elem/cl_DF_minusp.cc"

namespace cln {

CL_INLINE2 const cl_DF CL_INLINE2_DECL(float_sign) (const cl_DF& x)
{
// Methode: x>=0 -> Ergebnis 1.0; x<0 -> Ergebnis -1.0
	return (!minusp_inline(x) ? cl_DF_1 : cl_DF_minus1);
}

}  // namespace cln
