// integer_decode_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"

#include "base/cl_inline.h"
#include "float/sfloat/misc/cl_SF_idecode.cc"
#include "float/ffloat/misc/cl_FF_idecode.cc"
#include "float/dfloat/misc/cl_DF_idecode.cc"
#include "float/lfloat/misc/cl_LF_idecode.cc"

namespace cln {

const cl_idecoded_float CL_FLATTEN integer_decode_float (const cl_F& x)
{
	floatcase(x
	,	return integer_decode_float_inline(x);
	,	return integer_decode_float_inline(x);
	,	return integer_decode_float_inline(x);
	,	return integer_decode_float_inline(x);
	);
}

}  // namespace cln
