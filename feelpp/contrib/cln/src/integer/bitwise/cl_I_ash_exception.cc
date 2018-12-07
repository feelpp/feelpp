// ash_exception().

// General includes.
#include "base/cl_sysdep.h"

// Specification.
#include "cln/integer.h"


// Implementation.
#include "cln/io.h"
#include "cln/integer_io.h"
#include <sstream>

namespace cln {

static inline const std::string
ash_error_msg (const cl_I& badamount)
{
	std::ostringstream buf;
	fprint(buf, "ash: too large shift amount: ");
	fprint(buf, badamount);
	return buf.str();
}

ash_exception::ash_exception (const cl_I& badamount)
	: runtime_exception(ash_error_msg(badamount))
{}

}  // namespace cln
