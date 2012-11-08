// binary operator -

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "cln/rational.h"
#include "cln/float.h"

namespace cln {

const cl_R operator- (const cl_R& x, const cl_R& y)
{
	if (eq(y,0)) { return x; }
	elif (eq(x,0)) { return -y; }
	else
#define minus(a,b) a-b
GEN_R_OP2_2(x,y, minus, return)
}

}  // namespace cln
