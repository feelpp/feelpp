// catalanconst().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"
#include "float/transcendental/cl_F_tran.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F catalanconst (const cl_F& y)
{
	floattypecase(y
	,	return cl_SF_catalanconst();
	,	return cl_FF_catalanconst();
	,	return cl_DF_catalanconst();
	,	return catalanconst(TheLfloat(y)->len);
	);
}

}  // namespace cln
