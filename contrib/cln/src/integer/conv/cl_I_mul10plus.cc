// mul_10_plus_x().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

const cl_I mul_10_plus_x (const cl_I& y, unsigned char x)
{
	CL_ALLOCA_STACK;
	var uintD* MSDptr;
	var uintC len;
	var uintD* LSDptr;
	I_to_NDS_1(y, MSDptr=,len=,LSDptr=); // NDS zu Y
	var uintD carry = mulusmall_loop_lsp(10,LSDptr,len,x); // mal 10, plus x
	if (!(carry==0))
		{ lsprefnext(MSDptr) = carry; len++; }
	return UDS_to_I(MSDptr,len); // UDS als Integer zurück
}

}  // namespace cln
