// cl_ln10().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "float/transcendental/cl_F_tran.h"


// Implementation.

#include "cln/lfloat.h"
#include "float/lfloat/cl_LF.h"

namespace cln {

static inline const cl_LF compute_ln10_old (uintC len)
{
	return ln(cl_I_to_LF(10,len));
}

// ln 10 =
// = 46 atanh(1/31) + 34 atanh(1/49) + 20 atanh(1/161)
// = 478 atanh(1/251) + 180 atanh(1/449) - 126 atanh(1/4801) + 206 atanh(1/8749)

static inline const cl_LF compute_ln10_p235 (uintC len)
{
	var uintC actuallen = len+1;
	return shorten(  The(cl_LF)(46 * cl_atanh_recip(31,actuallen))
	               + The(cl_LF)(34 * cl_atanh_recip(49,actuallen))
	               + The(cl_LF)(20 * cl_atanh_recip(161,actuallen)),
	               len
	              );
}

static inline const cl_LF compute_ln10_p2357 (uintC len)
{
	var uintC actuallen = len+1;
	return shorten(  The(cl_LF)(478 * cl_atanh_recip(251,actuallen))
	               + The(cl_LF)(180 * cl_atanh_recip(449,actuallen))
	               - The(cl_LF)(126 * cl_atanh_recip(4801,actuallen))
	               + The(cl_LF)(206 * cl_atanh_recip(8749,actuallen)),
	               len
	              );
}

#define compute_ln10 compute_ln10_p2357

const cl_LF cl_ln10 (uintC len)
{
	var uintC oldlen = TheLfloat(cl_LF_ln10())->len; // vorhandene Länge
	if (len < oldlen)
		return shorten(cl_LF_ln10(),len);
	if (len == oldlen)
		return cl_LF_ln10();

	// TheLfloat(cl_LF_ln10())->len um mindestens einen konstanten Faktor
	// > 1 wachsen lassen, damit es nicht zu häufig nachberechnet wird:
	var uintC newlen = len;
	oldlen += floor(oldlen,2); // oldlen * 3/2
	if (newlen < oldlen)
		newlen = oldlen;

	// gewünschte > vorhandene Länge -> muß nachberechnen:
	cl_LF_ln10() = compute_ln10(newlen);
	return (len < newlen ? shorten(cl_LF_ln10(),len) : cl_LF_ln10());
}

}  // namespace cln
