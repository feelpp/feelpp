// float_sign().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/lfloat.h"


// Implementation.

#include "float/lfloat/cl_LF.h"
#include "float/lfloat/cl_LF_impl.h"

#include "base/cl_inline.h"
#include "float/lfloat/elem/cl_LF_minusp.cc"

namespace cln {

CL_INLINE2 const cl_LF CL_INLINE2_DECL(float_sign) (const cl_LF& x)
{
// Methode: x>=0 -> Ergebnis 1.0; x<0 -> Ergebnis -1.0
	return encode_LF1s(TheLfloat(x)->sign,TheLfloat(x)->len);
}

}  // namespace cln
