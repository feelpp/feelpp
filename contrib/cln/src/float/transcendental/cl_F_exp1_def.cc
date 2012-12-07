// exp1().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "float/cl_F.h"
#include "float/transcendental/cl_F_tran.h"

namespace cln {

const cl_F exp1 (void)
{
	floatformatcase(default_float_format
	,	return cl_SF_exp1();
	,	return cl_FF_exp1();
	,	return cl_DF_exp1();
	,	return exp1(len);
	);
}

}  // namespace cln
