// recip2adic().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "base/digitseq/cl_2DS.h"


// Implementation.

#include "base/digit/cl_2D.h"
#include "base/digitseq/cl_DS.h"

// Break-even point of the Newton iteration vs. standard div2adic.
#if CL_USE_GMP
const unsigned int recip2adic_threshold = 620;
#else
// Use the old default values from CLN version <= 1.0.3 as a crude estimate.
const unsigned int recip2adic_threshold = 380;
#endif

namespace cln {

void recip2adic (uintC len, const uintD* a_LSDptr, uintD* dest_LSDptr)
{
	// Method:
	// If len < threshold, use regular 2-adic division.
	// Else [Newton iteration] set n := ceiling(len/2),
	//   compute recursively b := recip2adic(a mod 2^(intDsize*n)),
	//   return 2*b-a*b^2 mod 2^(intDsize*2*n).
	CL_ALLOCA_STACK;
	var uintL k = 0; // number of Newton steps
	var uintC n = len;
	while (n >= recip2adic_threshold) {
		n = ceiling(n,2);
		k++;
	}
	// Nonrecursive step.
	var uintD* one_LSDptr;
	num_stack_alloc(n,,one_LSDptr=);
	lspref(one_LSDptr,0) = 1;
	clear_loop_lsp(one_LSDptr lspop 1,n-1);
	div2adic(n,one_LSDptr,a_LSDptr,dest_LSDptr);
	// Newton iteration.
	if (k > 0) {
		var uintD* b2_LSDptr;
		var uintD* prod_LSDptr;
		num_stack_alloc(len+1,,b2_LSDptr=);
		num_stack_alloc(2*len,,prod_LSDptr=);
		do {
			// n = ceiling(len/2^k)
			// Compute n2 = ceiling(len/2^(k-1)),
			// then n = ceiling(n2/2).
			k--;
			var uintC n2 = ((len-1)>>k)+1; // = 2*n or = 2*n-1
			// Set b := 2*b-a*b^2 mod 2^(intDsize*n2)
			cl_UDS_mul_square(dest_LSDptr,n,b2_LSDptr); // b^2
			cl_UDS_mul(b2_LSDptr,n2,a_LSDptr,n2,prod_LSDptr); // a*b^2
			clear_loop_lsp(dest_LSDptr lspop n,n2-n);
			shift1left_loop_lsp(dest_LSDptr,n+1); // (n+1 instead of n2 is ok)
			subfrom_loop_lsp(prod_LSDptr,dest_LSDptr,n2);
			n = n2;
		} while (k > 0);
	}
}
// Bit complexity (N := len): O(M(N)).

}  // namespace cln
