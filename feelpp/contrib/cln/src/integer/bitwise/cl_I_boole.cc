// boole().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.

namespace cln {

const cl_I boole (cl_boole op, const cl_I& x, const cl_I& y)
{
	switch (op) {
		case boole_clr:
			return 0;
		case boole_set:
			return -1;
		case boole_1:
			return x;
		case boole_2:
			return y;
		case boole_c1:
			return lognot(x);
		case boole_c2:
			return lognot(y);
		case boole_and:
			return logand(x,y);
		case boole_ior:
			return logior(x,y);
		case boole_xor:
			return logxor(x,y);
		case boole_eqv:
			return logeqv(x,y);
		case boole_nand:
			return lognand(x,y);
		case boole_nor:
			return lognor(x,y);
		case boole_andc1:
			return logandc1(x,y);
		case boole_andc2:
			return logandc2(x,y);
		case boole_orc1:
			return logorc1(x,y);
		case boole_orc2:
			return logorc2(x,y);
		default:
			NOTREACHED
	}
}

}  // namespace cln
