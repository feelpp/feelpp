// operator<< on cl_MI.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/modinteger.h"


// Implementation.

#include "cln/integer.h"

namespace cln {

const cl_MI operator<< (const cl_MI& x, sintC y) // assume 0 <= y < 2^(intCsize-1)
{
	if (y == 0)
		return x;
	if (y == 1) // frequent case
		return x+x;
	var const cl_modint_ring& R = x.ring();
	// Method:
	// Algorithm 1: divide (x.rep << y) by m.
	//              asymptotical cost: O(y * log m).
	// Algorithm 2: x * expt(2 mod m,y) using modular integer operations.
	//              asymptotical cost: O(log y * (log m)^2).
	// Use algorithm 1 for small y, algorithm 2 for large y.
	if ((R->bits < 0) || (y <= 2*R->bits))
		return cl_MI(R, R->reduce_modulo(x.rep << y));
	else
		return x * expt_pos(R->canonhom(2), (cl_I)(long)y);
}

}  // namespace cln
