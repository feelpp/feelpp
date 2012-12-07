// cl_trialdivision().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "numtheory/cl_IF.h"


// Implementation.

#include "integer/cl_I.h"

#if !(intDsize >= 16)
#error "intDsize too small for trialdivision!"
#endif

namespace cln {

uint32 cl_trialdivision (const cl_I& n, uint32 d1, uint32 d2)
{
	var uintL i = cl_small_prime_table_search(d1);
	var const uint16 * ptr = &cl_small_prime_table[i];
	var const uint16 * ptr_limit = &cl_small_prime_table[cl_small_prime_table_search(d2+1)];
	// Unpack n.
	CL_ALLOCA_STACK;
	var const uintD* n_MSDptr;
	var uintC n_len;
	I_to_NDS_nocopy(n, n_MSDptr=,n_len=,,false,);
	if (mspref(n_MSDptr,0)==0) { msshrink(n_MSDptr); n_len--; }
	// Make room for a quotient.
	var uintD* q_MSDptr;
	num_stack_alloc(n_len,q_MSDptr=,);
	// Division loop.
	for ( ; ptr < ptr_limit; ptr++) {
		var uint32 prime = *ptr;
		var uintD r = divucopy_loop_msp(prime,n_MSDptr,q_MSDptr,n_len);
		if (r == 0)
			return prime;
	}
	return 0;
}

}  // namespace cln
