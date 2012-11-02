// cl_I_LF_div().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/lfloat/cl_LF.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/lfloat/cl_LF_impl.h"
#include "cln/integer.h"
#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"
#include "float/cl_F.h"
#include "base/cl_N.h"

namespace cln {

const cl_R cl_I_LF_div (const cl_I& x, const cl_LF& y)
{
// Method:
// If x=0, return 0.
// Else convert x to a float and divide.
// (If x is shorter than y, we would gain nothing by dividing the absolute
// value of x by the mantissa of y, since the numerator of the division would
// have to have 2*length(y)+1 words, even if length(x) is much smaller than
// length(y).)
	if (eq(x,0)) { return 0; }
	var uintC len = TheLfloat(y)->len;
	return cl_I_to_LF(x,len) / y;
}
// Bit complexity (N = max(length(x),length(y))): O(M(N)).

}  // namespace cln
