// read_rational().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/rational_io.h"


// Implementation.

#include "rational/cl_RA.h"
#include "integer/cl_I.h"
#include "base/cl_N.h"

namespace cln {

const cl_RA read_rational (unsigned int base, cl_signean sign, const char * string, uintC index1, uintC index3, uintC index2)
{
	var uintC index3_1 = index3+1; // Index der ersten Nennerziffer
	var cl_I den = // Nenner
	  digits_to_I(&string[index3_1],index2-index3_1,(uintD)base);
	if (zerop(den)) // Division durch 0 abfangen
		{ throw division_by_0_exception(); }
	var cl_I num = // Zähler
	  digits_to_I(&string[index1],index3-index1,(uintD)base);
	if (!(sign == 0))
		num = -num; // incl. Vorzeichen
	return I_posI_div_RA(num,den);
}

}  // namespace cln
