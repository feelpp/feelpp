// cl_ln2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "float/cl_F.h"

namespace cln {

const cl_F cl_ln2 (float_format_t f)
{
	floatformatcase((uintC)f
	,	return cl_SF_ln2();
	,	return cl_FF_ln2();
	,	return cl_DF_ln2();
	,	return cl_ln2(len);
	);
}

}  // namespace cln
