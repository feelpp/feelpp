// read_integer().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer_io.h"


// Implementation.

#include "integer/cl_I.h"

namespace cln {

const cl_I read_integer (unsigned int base, cl_signean sign, const char * string, uintC index1, uintC index2)
{
	var cl_I x = digits_to_I(&string[index1],index2-index1,(uintD)base);
	if (sign == 0)
		return x;
	else
		return -x; // negatives Vorzeichen -> Vorzeichenwechsel
}

}  // namespace cln
