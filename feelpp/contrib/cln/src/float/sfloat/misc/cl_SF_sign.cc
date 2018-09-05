// float_sign().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_minusp.cc"

namespace cln {

CL_INLINE2 const cl_SF CL_INLINE2_DECL(float_sign) (const cl_SF& x)
{
// Methode: x>=0 -> Ergebnis 1.0; x<0 -> Ergebnis -1.0
	return encode_SF(SF_sign(x),1,bit(SF_mant_len));
}

}  // namespace cln
