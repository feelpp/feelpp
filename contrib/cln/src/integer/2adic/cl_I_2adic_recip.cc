// cl_recip2adic().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "base/digitseq/cl_DS.h"
#include "base/digitseq/cl_2DS.h"
#include "integer/bitwise/cl_I_log.h"

namespace cln {

const cl_I cl_recip2adic (uintL n, const cl_I& x)
{
	var uintL len = ceiling(n,intDsize);
	CL_ALLOCA_STACK;
	var const uintD* x_LSDptr;
	if (bignump(x) && TheBignum(x)->length >= len)
		// no need to copy x
		x_LSDptr = BN_LSDptr(x);
	else {	// copy x
		var uintL x_len = I_to_DS_need(x);
		if (x_len < len) { x_len = len; }
		I_to_DS_n(x,x_len,x_LSDptr=);
		x_LSDptr = x_LSDptr mspop x_len;
	}
	var uintD* y_LSDptr;
	num_stack_alloc_1(len,,y_LSDptr=);
	// Compute inverse mod 2^(intDsize*len).
	recip2adic(len,x_LSDptr,y_LSDptr);
	// Reduce mod 2^n.
	if ((n % intDsize) != 0)
		lspref(y_LSDptr,floor(n,intDsize)) &= (bit(n % intDsize) - 1);
	return UDS_to_I(y_LSDptr lspop len,len);
}
// Bit complexity (N := n): O(M(N)).

}  // namespace cln
