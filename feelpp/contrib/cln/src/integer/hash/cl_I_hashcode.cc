// cl_I hashcode().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "integer/cl_I.h"


// Implementation.

namespace cln {

uintptr_t hashcode (const cl_I& x)
{
	var uintptr_t code = 0x814BE3A5;
	// We walk through all limbs. It may take some time for very large
	// integers, but it's better than completely ignoring some limbs.
	if (fixnump(x)) {
		#if (cl_value_len <= intLsize)
		code += FN_to_V(x);
		#elif (cl_word_size==64)
		code += FN_to_Q(x);
		code ^= (code >> 32);
		#endif
		code &= 0xFFFFFFFF;
	} else {
		var const uintD* MSDptr;
		var uintC len;
		BN_to_NDS_nocopy(x, MSDptr=,len=,);
		for (; len > 0; len--) {
			var uintD c = msprefnext(MSDptr);
			code = (code << 5) | (code >> 27); // rotate left 5 bits
			code += (intptr_t)c << 16;
			code ^= (intptr_t)c;
			code &= 0xFFFFFFFF;
		}
	}
	return code;
}

}  // namespace cln
