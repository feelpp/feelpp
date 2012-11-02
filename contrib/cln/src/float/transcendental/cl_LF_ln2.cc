// cl_ln2().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

static inline const cl_LF compute_ln2_old (uintC len)
{
	// Here, it is tricky to avoid a recursive loop. We assume that ln()
	// will not invoke cl_ln2() if its argument is between 2/3 and 4/3.
	// So we compute -2*ln(1/sqrt(2)).
	return -scale_float(ln(sqrt(scale_float(cl_I_to_LF(1,len),-1))),1);
}

// ln 2 =
//  = 2 atanh(1/3)
//  = 4 atanh(1/7) + 2 atanh(1/17)
//  = 14 atanh(1/31) + 10 atanh(1/49) + 6 atanh(1/161)
//  = 144 atanh(1/251) + 54 atanh(1/449) - 38 atanh(1/4801) + 62 atanh(1/8749)

static inline const cl_LF compute_ln2_p2 (uintC len)
{
	return scale_float(cl_atanh_recip(3,len),1);
}

static inline const cl_LF compute_ln2_p23 (uintC len)
{
	var uintC actuallen = len+1;
	return shorten(scale_float(cl_atanh_recip(7,actuallen),2)
	               + scale_float(cl_atanh_recip(17,actuallen),1),
	               len
	              );
}

static inline const cl_LF compute_ln2_p235 (uintC len)
{
	var uintC actuallen = len+1;
	return shorten(  The(cl_LF)(14 * cl_atanh_recip(31,actuallen))
	               + The(cl_LF)(10 * cl_atanh_recip(49,actuallen))
	               + The(cl_LF)( 6 * cl_atanh_recip(161,actuallen)),
	               len
	              );
}

static inline const cl_LF compute_ln2_p2357 (uintC len)
{
	var uintC actuallen = len+1;
	return shorten(  The(cl_LF)(144 * cl_atanh_recip(251,actuallen))
	               + The(cl_LF)( 54 * cl_atanh_recip(449,actuallen))
	               - The(cl_LF)( 38 * cl_atanh_recip(4801,actuallen))
	               + The(cl_LF)( 62 * cl_atanh_recip(8749,actuallen)),
	               len
	              );
}

// Timings of the above algorithms, on an i486 33 MHz, running Linux.
// N = 250
// p2       1.71 s                        total 1.71 s
// p23      0.88 + 0.60 s                 total 1.48 s
// p235     0.51 + 0.48 + 0.40 s          total 1.39 s
// p2357    0.38 + 0.37 + 0.31 + 0.29 s   total 1.35 s
// In general, for fixed precision N, the computation time of atanh(1/m)
// seems to be proportional to 1/log(m).

#define compute_ln2 compute_ln2_p2357

const cl_LF cl_ln2 (uintC len)
{
	var uintC oldlen = TheLfloat(cl_LF_ln2())->len; // vorhandene Länge
	if (len < oldlen)
		return shorten(cl_LF_ln2(),len);
	if (len == oldlen)
		return cl_LF_ln2();

	// TheLfloat(cl_LF_ln2())->len um mindestens einen konstanten Faktor
	// > 1 wachsen lassen, damit es nicht zu häufig nachberechnet wird:
	var uintC newlen = len;
	oldlen += floor(oldlen,2); // oldlen * 3/2
	if (newlen < oldlen)
		newlen = oldlen;

	// gewünschte > vorhandene Länge -> muß nachberechnen:
	cl_LF_ln2() = compute_ln2(newlen);
	return (len < newlen ? shorten(cl_LF_ln2(),len) : cl_LF_ln2());
}

}  // namespace cln
