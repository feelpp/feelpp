// cl_F_extendsqrtx().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/cl_F.h"


// Implementation.

#include "float/sfloat/cl_SF.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_F cl_F_extendsqrtx (const cl_F& x)
{
// Methode:
// SF -> DF wegen 17+sqrt(17)+2+7 = 30.2 < 53
// FF -> DF wegen 24+sqrt(24)+2+7 = 37.9 < 53
// DF -> LF(5) wegen 53+sqrt(53)+2+10 = 72.3 < 80
// LF(n) -> LF(n+i)
	floatcase(x
	,	return cl_SF_to_DF(x); // 17+sqrt(17)+2+7 = 30.2 < 53
	,	return cl_FF_to_DF(x); // 24+sqrt(24)+2+7 = 37.9 < 53
	,	return cl_DF_to_LF(x,ceiling(73,intDsize)); // 53+sqrt(53)+2+10 = 72.3 < 73
	,	return extend(x,cl_LF_len_incsqrtx(TheLfloat(x)->len));
	);
}

}  // namespace cln
