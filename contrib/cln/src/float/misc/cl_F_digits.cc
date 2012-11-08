// float_digits().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/misc/cl_SF_digits.cc"
#include "float/ffloat/misc/cl_FF_digits.cc"
#include "float/dfloat/misc/cl_DF_digits.cc"
#include "float/lfloat/misc/cl_LF_digits.cc"

namespace cln {

uintC float_digits (const cl_F& x)
{
	floatcase(x
	,	return float_digits_inline(x);
	,	return float_digits_inline(x);
	,	return float_digits_inline(x);
	,	return float_digits_inline(x);
	);
}

}  // namespace cln
