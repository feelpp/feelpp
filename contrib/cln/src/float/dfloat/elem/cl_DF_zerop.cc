// zerop().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

CL_INLINE bool CL_INLINE_DECL(zerop) (const cl_DF& x)
{
#if 0
	return DF_uexp(TheDfloat(x)->dfloat_value_semhi) == 0;
#else // this is faster
	return TheDfloat(x)->dfloat_value_semhi == 0;
#endif
}

}  // namespace cln
