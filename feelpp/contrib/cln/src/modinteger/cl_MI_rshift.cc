// operator>> on cl_MI.

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/modinteger.h"


// Implementation.

#include "cln/integer.h"
#include "cln/exception.h"
#include "modinteger/cl_MI.h"

namespace cln {

const cl_MI operator>> (const cl_MI& x, sintC y) // assume 0 <= y < 2^(intCsize-1)
{
	if (y == 0)
		return x;
	var const cl_modint_ring& R = x.ring();
	if (!oddp(R->modulus)) {
		if (R->modulus == 2)
			throw division_by_0_exception();
		else
			return (cl_MI_x)cl_notify_composite(R,2);
	}
	if (y == 1) // frequent case
		return cl_MI(R, (evenp(x.rep) ? x.rep : x.rep + R->modulus) >> 1);
	// Method:
	// Algorithm 1: add a multiple of m to x.rep so that it becomes
	//              divisible by 2^y (2-adic division), then shift right.
	//              asymptotical cost: O(y * log m).
	// Algorithm 2: x * expt(2 mod m,-y) using modular integer operations.
	//              asymptotical cost: O(log y * (log m)^2).
	// Use algorithm 1 for small y, algorithm 2 for large y.
#if 0
	if (y <= 2*R->bits)
		throw runtime_exception(); // not yet implemented
	else
#endif
		return R->div(x, expt_pos(R->canonhom(2), (cl_I)(long)y));
}

}  // namespace cln
