// NUDS_to_I().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

#include "cln/number.h"
#include "base/digitseq/cl_DS.h"

#include "base/cl_inline.h"
#include "integer/conv/cl_I_from_NDS.cc"

namespace cln {

CL_INLINE2 const cl_I CL_INLINE2_DECL(NUDS_to_I) (uintD* MSDptr, uintC len)
{
	if ((!(len==0)) && ((sintD)mspref(MSDptr,0) < 0))
		// Falls die Länge >0 und das Most significant Bit = 1 sind,
		// die Digit Sequence um ein Nulldigit erweitern:
		{ lsprefnext(MSDptr) = 0; len++; }
	return NDS_to_I_inline(MSDptr,len);
}

}  // namespace cln
