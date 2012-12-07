// cl_div2adic().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "base/digitseq/cl_DS.h"
#include "base/digitseq/cl_2DS.h"
#include "integer/bitwise/cl_I_log.h"

namespace cln {

const cl_I cl_div2adic (uintL n, const cl_I& x, const cl_I& y)
{
	var uintL len = ceiling(n,intDsize);
	CL_ALLOCA_STACK;
	var const uintD* x_LSDptr;
	var const uintD* y_LSDptr;
	if (bignump(x) && TheBignum(x)->length >= len)
		// no need to copy x
		x_LSDptr = BN_LSDptr(x);
	else {	// copy x
		var uintL x_len = I_to_DS_need(x);
		if (x_len < len) { x_len = len; }
		I_to_DS_n(x,x_len,x_LSDptr=);
		x_LSDptr = x_LSDptr mspop x_len;
	}
	if (bignump(y) && TheBignum(y)->length >= len)
		// no need to copy y
		y_LSDptr = BN_LSDptr(y);
	else {	// copy y
		var uintL y_len = I_to_DS_need(y);
		if (y_len < len) { y_len = len; }
		I_to_DS_n(y,y_len,y_LSDptr=);
		y_LSDptr = y_LSDptr mspop y_len;
	}
	var uintD* z_LSDptr;
	num_stack_alloc_1(len,,z_LSDptr=);
	// Compute quotient mod 2^(intDsize*len).
	div2adic(len,x_LSDptr,y_LSDptr,z_LSDptr);
	// Reduce mod 2^n.
	if ((n % intDsize) != 0)
		lspref(z_LSDptr,floor(n,intDsize)) &= (bit(n % intDsize) - 1);
	return UDS_to_I(z_LSDptr lspop len,len);
}
// Bit complexity (N := n): O(M(N)).

}  // namespace cln
