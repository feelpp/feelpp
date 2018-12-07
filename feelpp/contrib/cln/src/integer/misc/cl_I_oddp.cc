// oddp().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

#include "integer/cl_I.h"
#include "base/digitseq/cl_DS.h"

namespace cln {

bool oddp (const cl_I& x)
{
	if (fixnump(x)) {
		// Fixnum: Bit 0 abprüfen
		if (x.word & bit(cl_value_shift))
			return true;
		else
			return false;
	} else {
		// Bignum: Bit 0 im letzten Digit abprüfen
		if (lspref(BN_LSDptr(x),0) & bit(0))
			return true;
		else
			return false;
	}
}

}  // namespace cln
