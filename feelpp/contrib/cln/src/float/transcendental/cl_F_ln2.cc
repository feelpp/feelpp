// cl_ln2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "float/cl_F.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F cl_ln2 (const cl_F& y)
{
	floattypecase(y
	,	return cl_SF_ln2();
	,	return cl_FF_ln2();
	,	return cl_DF_ln2();
	,	return cl_ln2(TheLfloat(y)->len);
	);
}

}  // namespace cln
