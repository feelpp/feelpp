// unary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/sfloat.h"


// Implementation.

#include "float/sfloat/cl_SF.h"

#include "base/cl_inline.h"
#include "float/sfloat/elem/cl_SF_zerop.cc"

namespace cln {

const cl_SF operator- (const cl_SF& x)
{
// Methode:
// Falls x=0.0, fertig. Sonst Vorzeichenbit umdrehen.
	if (zerop_inline(x))
		return SF_0;
	return cl_SF_from_word(x.word ^ ((cl_uint)1 << SF_sign_shift));
}

}  // namespace cln
