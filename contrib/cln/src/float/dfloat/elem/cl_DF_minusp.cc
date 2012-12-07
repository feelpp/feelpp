// minusp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/dfloat.h"


// Implementation.

#include "float/dfloat/cl_DF.h"

namespace cln {

CL_INLINE bool CL_INLINE_DECL(minusp) (const cl_DF& x)
{
#if (cl_word_size==64)
	return (sint64)TheDfloat(x)->dfloat_value_semhi < 0;
#else
	return (sint32)TheDfloat(x)->dfloat_value_semhi < 0;
#endif
}

}  // namespace cln
