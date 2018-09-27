// random_R().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/real.h"


// Implementation.

#include "real/cl_R.h"
#include "rational/cl_RA.h"
#include "cln/io.h"
#include "cln/real_io.h"
#include "cln/exception.h"
#include <sstream>

namespace cln {

const cl_R random_R (random_state& r, const cl_R& n)
{
	// n muß eine reelle Zahl sein, >0 und Float oder Integer
	if (plusp(n)) {
		if (floatp(n)) {
			DeclareType(cl_F,n);
			return random_F(r,n);
		} else {
			DeclareType(cl_RA,n);
			if (integerp(n)) {
				DeclareType(cl_I,n);
				return random_I(r,n);
			}
		}
	}
	std::ostringstream buf;
	fprint(buf, "random: argument should be positive and an integer or float: ");
	fprint(buf, n);
	throw runtime_exception(buf.str());
}

}  // namespace cln
