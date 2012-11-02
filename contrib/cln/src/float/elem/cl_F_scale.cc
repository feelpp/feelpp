// scale_float().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/float.h"


// Implementation.

#include "cln/sfloat.h"
#include "cln/ffloat.h"
#include "cln/dfloat.h"
#include "cln/lfloat.h"
#include "float/cl_F.h"

namespace cln {

const cl_F scale_float (const cl_F& x, sintC delta)
{
	floatcase(x
	,	return scale_float(x,delta);
	,	return scale_float(x,delta);
	,	return scale_float(x,delta);
	,	return scale_float(x,delta);
	);
}

}  // namespace cln
