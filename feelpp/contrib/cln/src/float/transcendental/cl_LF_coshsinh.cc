// cl_coshsinh_ratseries().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"
#include "cln/integer.h"

namespace cln {

inline const cl_LF_cosh_sinh_t operator* (const cl_LF_cosh_sinh_t& a, const cl_LF_cosh_sinh_t& b)
{
	return cl_LF_cosh_sinh_t(a.cosh*b.cosh+a.sinh*b.sinh,a.sinh*b.cosh+a.cosh*b.sinh);
}

const cl_LF_cosh_sinh_t cl_coshsinh_ratseries (const cl_LF& x)
{
	// Similar to expx_ratseries.
	var uintC len = TheLfloat(x)->len;
	var cl_idecoded_float x_ = integer_decode_float(x);
	// x = (-1)^sign * 2^exponent * mantissa
	var uintE lq = cl_I_to_UE(- x_.exponent);
	var const cl_I& p = x_.mantissa;
	// Compute sinh(p/2^lq) and cosh(p/2^lq) by splitting into pieces.
	var bool first_factor = true;
	var cl_LF_cosh_sinh_t product;
	var uintE b1;
	var uintE b2;
	for (b1 = 0, b2 = 1; b1 < lq; b1 = b2, b2 = 2*b2) {
		// Piece containing bits b1+1..b2 after "decimal point"
		// in the binary representation of (p/2^lq).
		var uintE lqk = (lq >= b2 ? b2 : lq);
		var cl_I pk = ldb(p,cl_byte(lqk-b1,lq-lqk));
		// Compute sinh(pk/2^lqk) and cosh(pk/2^lqk).
		if (!zerop(pk)) {
			if (minusp(x_.sign)) { pk = -pk; }
			var cl_LF_cosh_sinh_t factor = cl_coshsinh_aux(pk,lqk,len);
			if (first_factor) {
				product = factor;
				first_factor = false;
			} else
				product = product * factor;
		}
	}
	if (first_factor)
		return cl_LF_cosh_sinh_t(cl_I_to_LF(1,len),cl_I_to_LF(0,len));
	else
		return product;
}
// Bit complexity (N = length(x)): O(log(N)^2*M(N)).

}  // namespace cln
