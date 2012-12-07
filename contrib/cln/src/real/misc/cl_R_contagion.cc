// contagion().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "real/cl_R.h"


// Implementation.

#include "float/lfloat/cl_LF.h"

namespace cln {

const cl_R contagion (const cl_R& x, const cl_R& y)
{
#define X  { return x; }
#define Y  { return y; }
#if 0
	if (R_rationalp(x)) Y
	elif (R_rationalp(y)) X
	else
		floattypecase(x
		, X // floattypecase(y, X,X,X,X)
		, floattypecase(y, Y,X,X,X)
		, floattypecase(y, Y,Y,X,X)
		, floattypecase(y, Y,Y,Y,
			if (TheLfloat(x)->len <= TheLfloat(y)->len) X else Y
			)
		);
#else // faster type dispatch
	realtypecase(x
	,	Y
	,	Y
	,	Y
	,	X // realtypecase(y, X,X,X, X,X,X,X)
	,	realtypecase(y, X,X,X, Y,X,X,X)
	,	realtypecase(y, X,X,X, Y,Y,X,X)
	,	realtypecase(y, X,X,X, Y,Y,Y,
			if (TheLfloat(x)->len <= TheLfloat(y)->len) X else Y
			)
	);
#endif
}

}  // namespace cln
